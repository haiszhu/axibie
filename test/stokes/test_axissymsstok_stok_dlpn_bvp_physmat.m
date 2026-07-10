clearvars; format short e;
addpath('/Users/hzhu/Documents/Github/axibie/utils');
addpath('/Users/hzhu/Documents/Github/AxiStokes3D/utils');
addpath('/Users/hzhu/Documents/Github/AxiStokes3D/matlab');
addpath('/Users/hzhu/Documents/Github/AxiStokes3D/external/kdtree/toolbox')
addpath('/Users/hzhu/Documents/Github/AxiStokes3D/external/fmm3d/matlab');   % stfmm3d/Sto3d{SLP,DLP}fmm

% ---- geometry switch (mirror the SLPn physmat2 test) ----
shape = 'sphere';                                   % 'sphere' | 'ellipse' | 'cshape'
p=16; iside=1; iclosed=0; mu=1; gate=2.0; fmm_eps=1e-14;
if strcmp(shape,'sphere')
  Z=@(t) sin(t)-1i*cos(t);                          % unit sphere
  np_vals=2:16;
  y1=[0.30;0;0.20]; y2=[0.32*cos(pi/4);0.32*sin(pi/4);-0.25];
  gv=linspace(-1.9,1.9,32);
elseif strcmp(shape,'ellipse')
  a=2; Z=@(t) a^(-1/3)*sin(t)-1i*a^(2/3)*cos(t);    % volume-conserving prolate ellipsoid
  np_vals=2:16;
  y1=[0.20;0;0.40]; y2=[0.25*cos(pi/4);0.25*sin(pi/4);-0.50];
  gv=linspace(-2.4,2.4,18);
elseif strcmp(shape,'cshape')
  lam=0.75; Z=@(t) -(1.5+cos(t)).*(-sin(lam*pi*sin(t))+1i*cos(lam*pi*sin(t)));
  np_vals=6:2:24;
  y1=[0.70;0;-1.20]; y2=[1.30*cos(pi/4);1.30*sin(pi/4);0.20];
  gv=linspace(-3,3,24);
end
F1=[1;-0.7;0.5]; F2=[0.4;0.9;-0.6];                 % interior Stokeslet forces

% interior Stokeslets -> exact exterior Stokes velocity u_ex; Neumann data = traction on surface.
uex=@(X) stk(X,y1,F1)+stk(X,y2,F2);

% ---- exterior 3D target grid (once): outside the revolved meridian, >1e-2 off the surface ----
mer=Z(linspace(0,pi,400).');
[GX,GY,GZ]=meshgrid(gv,gv,gv); rho3=hypot(GX(:),GY(:)); z3=GZ(:); zt=rho3+1i*z3;
dmer=inf(size(zt)); for j=1:numel(mer), dmer=min(dmer,abs(zt-mer(j))); end
ext3=(~inpolygon(rho3,z3,real(mer),imag(mer))) & (dmer>1e-2);
P=[GX(ext3).';GY(ext3).';GZ(ext3).']; M3=size(P,2);
rho3=hypot(P(1,:),P(2,:)); z3=P(3,:); d3=dmer(ext3).';

errmax=nan(1,numel(np_vals));
err_trac=nan(1,numel(np_vals));

for kk=1:numel(np_vals)
  np=np_vals(kk); pmodes=2*np; nmodes=2*pmodes+1; nang=nmodes;

  % 1. geometry + Neumann data setup (LAB traction at ALL ring nodes, INTERLOCKED)
  s=[]; s.p=p; s.Z=Z; s.tpan=linspace(0,pi,np+1)'; s=quadr(s,[],'p','G');
  N=numel(s.x);
  sx=s.x(:); snx=s.nx(:); sws=s.ws(:); swxp=s.wxp(:); tpan=s.tpan(:);
  s3d=axisym_to_3d_quadrature(real(s.x),imag(s.x),s.ws,nang);
  phis=2*pi*(0:nang-1)/nang;
  s3dnx=[reshape(real(snx)*cos(phis),1,[]); reshape(real(snx)*sin(phis),1,[]); repmat(imag(snx).',1,nang)];
  nsf=N*nang;
  tlab=zeros(3,nsf);
  nr=real(snx).'; nz=imag(snx).';
  for i=1:nang
    th=2*pi*(i-1)/nang; X=s3d.x(:,(i-1)*N+(1:N)); N3=[nr*cos(th); nr*sin(th); nz];
    tlab(:,(i-1)*N+(1:N)) = stktrac(X,N3,y1,F1)+stktrac(X,N3,y2,F2);
  end
  tphys=tlab(:);

  % 2. loop over patch and compute close interaction nodes (self: azimuth-0 orbit reps)
  % FROZEN close-correction path (uncomment to compare): tx = sx; nt = N;
  % FROZEN close-correction path (uncomment to compare): t3dx  = [real(sx).'; zeros(1,N); imag(sx).'];
  % FROZEN close-correction path (uncomment to compare): tn3dx = [real(snx).'; zeros(1,N); imag(snx).'];      % on-surface orbit reps have n_theta=0
  % FROZEN close-correction path (uncomment to compare): tcxi = zeros(np+1,1); ntcx = 0;
  % FROZEN close-correction path (uncomment to compare): [tcxi, ntcx] = axps_closesize_mex(nt, tx, t3dx, p, np, sx, sws, gate, tcxi, ntcx);

  % 3. allocate space for the S' and D' physical-space close-correction matrices
  % FROZEN close-correction path (uncomment to compare): Ss = zeros(3*ntcx, 3*nang*p); ids = zeros(ntcx,1);
  % FROZEN close-correction path (uncomment to compare): Dd = zeros(3*ntcx, 3*nang*p); idd = zeros(ntcx,1);

  % 4. actual computation (S' + D' traction close corrections)
  % FROZEN close-correction path (uncomment to compare): [Ss, ids] = axps_closestokslpn_mex(nt, tx, t3dx, tn3dx, p, np, nang, sx, snx, sws, swxp, tpan, gate, ...
  % FROZEN close-correction path (uncomment to compare):     s3d.x, s3dnx, s3d.w, pmodes, iside, iclosed, mu, ntcx, tcxi, Ss, ids);
  % FROZEN close-correction path (uncomment to compare): [Dd, idd] = axps_closestokdlpn_mex(nt, tx, t3dx, tn3dx, p, np, nang, sx, snx, sws, swxp, tpan, gate, ...
  % FROZEN close-correction path (uncomment to compare):     s3d.x, s3dnx, s3d.w, pmodes, iside, iclosed, mu, ntcx, tcxi, Dd, idd);

  % 5. matrix-free operator: stfmm3d naive S' + D' at the sources, patched by sparse close corrections
  % FROZEN close-correction path (uncomment to compare): applyA = @(x) apply_tracfmm_close(x, N, p, np, nang, ntcx, tcxi, ids, Ss, idd, Dd, ...
  % FROZEN close-correction path (uncomment to compare):                                   s3d.x, s3dnx, s3d.w, mu, fmm_eps);

  % 5'. LEVEL-2 PHYSICAL master (axp_physmat_setup_mex): assembled dense SELF operators, LAB
  %     node-interleaved (== axp_physmat_mex's xA); combined -1/2 I + S' + D' (jump rides in S')
  nblk = 3*N*nang;
  A = axp_physmat_setup_mex(2,2,mu,1,1,p,np,pmodes,iside,iclosed,[1;N+1],[1;np+2],N,np+1, ...
      sx,snx,sws,swxp,tpan,0,[1;1],zeros(3,0),zeros(3,0),eye(3),zeros(3,1),nblk,nblk);       % -1/2 I + S'
  A = A + axp_physmat_setup_mex(2,4,mu,1,1,p,np,pmodes,iside,iclosed,[1;N+1],[1;np+2],N,np+1, ...
      sx,snx,sws,swxp,tpan,0,[1;1],zeros(3,0),zeros(3,0),eye(3),zeros(3,1),nblk,nblk);       % + D'

  % 6. solve the physical-space combined Neumann system (dense direct backslash, as the multi
  %    dense reference: the S' regularizes the hypersingular D' -> full-rank)
  % FROZEN close-correction path (uncomment to compare): nh = s3dnx(:)/norm(s3dnx(:));
  % FROZEN close-correction path (uncomment to compare): [sigphys,flag,relres,iter] = gmres(@(x) applyA(x), tphys, [], 1e-12, 200);
  % FROZEN close-correction path (uncomment to compare): fprintf('  gmres: flag=%d  iters=%d  relres=%.3e  sig.n=%.2e  ntcx=%d\n', flag, iter(2), relres, nh'*sigphys, ntcx);
  sigphys = A\tphys;

  % 7. kdtree near-target detection ONCE (ball per panel) -- the SAME near set feeds the
  %    velocity and traction verifications; far grid targets get pure FMM, no correction.
  axk_kdtree_build_mex(1, [rho3(:), z3(:), zeros(M3,1)], 64);
  isnear = false(M3,1);
  for k = 1:np
    jj = (k-1)*p + (1:p);
    qr = mean(real(s.x(jj))); qz = mean(imag(s.x(jj)));
    qradii = gate*sqrt(sum(s.ws(jj)));
    idxcin = axk_kdtree_ball_mex([qr, qz, 0], qradii, M3);
    isnear(idxcin) = true;
  end
  Pn = P(:,isnear); Mn = nnz(isnear);

  % 8. close sizes ONCE at the near set + velocity and traction compute passes
  rng(0); Nn = randn(3,Mn); Nn = Nn./vecnorm(Nn);
  ttn = (rho3(isnear) + 1i*z3(isnear)).';
  % FROZEN close-correction path (uncomment to compare): tcxi3 = zeros(np+1,1); ntcx3 = 0;
  % FROZEN close-correction path (uncomment to compare): [tcxi3, ntcx3] = axps_closesize_mex(Mn, ttn, Pn, p, np, sx, sws, gate, tcxi3, ntcx3);
  % FROZEN close-correction path (uncomment to compare): S3 = zeros(3*ntcx3, 3*nang*p); idx3 = zeros(ntcx3,1);
  % FROZEN close-correction path (uncomment to compare): [S3, idx3] = axps_closestokslp_mex(Mn, ttn, Pn, p, np, nang, sx, snx, sws, swxp, tpan, gate, ...
  % FROZEN close-correction path (uncomment to compare):     s3d.x, s3dnx, s3d.w, pmodes, iside, iclosed, mu, ntcx3, tcxi3, S3, idx3);
  % FROZEN close-correction path (uncomment to compare): D3 = zeros(3*ntcx3, 3*nang*p); id3 = zeros(ntcx3,1);
  % FROZEN close-correction path (uncomment to compare): [D3, id3] = closestokdlp_physmat(Mn, ttn, Pn, p, np, nang, sx, snx, sws, swxp, tpan, gate, ...
  % FROZEN close-correction path (uncomment to compare):     s3d.x, s3dnx, s3d.w, pmodes, iside, iclosed, mu, ntcx3, tcxi3, D3, id3);
  % FROZEN close-correction path (uncomment to compare): Stn = zeros(3*ntcx3, 3*nang*p); itn = zeros(ntcx3,1);
  % FROZEN close-correction path (uncomment to compare): [Stn, itn] = axps_closestokslpn_mex(Mn, ttn, Pn, Nn, p, np, nang, sx, snx, sws, swxp, tpan, gate, ...
  % FROZEN close-correction path (uncomment to compare):     s3d.x, s3dnx, s3d.w, pmodes, iside, iclosed, mu, ntcx3, tcxi3, Stn, itn);
  % FROZEN close-correction path (uncomment to compare): Dtn = zeros(3*ntcx3, 3*nang*p); idn = zeros(ntcx3,1);
  % FROZEN close-correction path (uncomment to compare): [Dtn, idn] = axps_closestokdlpn_mex(Mn, ttn, Pn, Nn, p, np, nang, sx, snx, sws, swxp, tpan, gate, ...
  % FROZEN close-correction path (uncomment to compare):     s3d.x, s3dnx, s3d.w, pmodes, iside, iclosed, mu, ntcx3, tcxi3, Dtn, idn);

  % 9'. matrix-free velocity eval: Sto3dSLPfmm + Sto3dDLPfmm at ALL targets + LEVEL-2 master eval
  %     (iinter=3, S + D) at the kdtree near set: the master fills its internal near zone (far rows
  %     stay 0), those rows replace the FMM values (naive full + near replacement)
  sigblk = reshape(reshape(sigphys,3,[]).',[],1);       % interlocked -> component-block
  ublk = Sto3dSLPfmm(struct('x',P), struct('x',s3d.x,'w',s3d.w), sigblk, fmm_eps) + ...
         Sto3dDLPfmm(struct('x',P), struct('x',s3d.x,'nx',s3dnx,'w',s3d.w), sigblk, fmm_eps);
  uvec = reshape(reshape(ublk,[],3).',[],1);            % component-block -> interlocked
  rown = reshape(3*(find(isnear).'-1)+(1:3)',[],1);
  un = uvec(rown);
  % FROZEN close-correction path (uncomment to compare): un = axps_closecorrapply_mex(3, Mn, p, np, nang, ntcx3, tcxi3, idx3, S3, sigphys, 0, un);
  % FROZEN close-correction path (uncomment to compare): un = axps_closecorrapply_mex(3, Mn, p, np, nang, ntcx3, tcxi3, id3, D3, sigphys, 0, un);
  Aev = axp_physmat_setup_mex(2,1,mu,3,1,p,np,pmodes,iside,iclosed,[1;N+1],[1;np+2],N,np+1, ...
      sx,snx,sws,swxp,tpan,Mn,[1;Mn+1],Pn,zeros(3,Mn),eye(3),zeros(3,1),3*Mn,nblk) ...
      + axp_physmat_setup_mex(2,3,mu,3,1,p,np,pmodes,iside,iclosed,[1;N+1],[1;np+2],N,np+1, ...
      sx,snx,sws,swxp,tpan,Mn,[1;Mn+1],Pn,zeros(3,Mn),eye(3),zeros(3,1),3*Mn,nblk);
  fillm = repelem(any(reshape(any(Aev~=0,2),3,Mn),1).',3,1);   % targets the master filled (interlocked rows)
  uev = Aev*sigphys; un(fillm) = uev(fillm);
  uvec(rown) = un;
  U3 = reshape(uvec, 3, M3);

  % 10. log10 velocity error (near = kdtree set, corrected; far = pure FMM)
  Uex3=uex(P); err3=vecnorm(U3-Uex3)/max(vecnorm(Uex3));
  errmax(kk)=max(err3(isnear));
  fprintf('Stokes DLPn+SLPn physical Neumann BVP [%s] (physmat, matrix-free): N=%d, np=%d, pmodes=%d, M3=%d\n',shape,N,np,pmodes,M3);
  fprintf('  kdtree-near (corrected): max err = %.3e  (%d pts)\n', max(err3(isnear)), Mn);
  fprintf('  far (pure FMM):          max err = %.3e  (%d pts)\n', max(err3(~isnear)), nnz(~isnear));

  % 11. S' + D' traction verification at the SAME near set (naive FMM base + master eval replacement)
  Tex = stktrac(Pn,Nn,y1,F1) + stktrac(Pn,Nn,y2,F2);
  tnear = slpntracfmm(sigphys, s3d.x, s3d.w, Pn, Nn, fmm_eps, 0) + ...
          dlpntracfmm(sigphys, s3d.x, s3dnx, s3d.w, Pn, Nn, mu, fmm_eps, 0);
  % FROZEN close-correction path (uncomment to compare): tnear = axps_closecorrapply_mex(3, Mn, p, np, nang, ntcx3, tcxi3, itn, Stn, sigphys, 0, tnear);
  % FROZEN close-correction path (uncomment to compare): tnear = axps_closecorrapply_mex(3, Mn, p, np, nang, ntcx3, tcxi3, idn, Dtn, sigphys, 0, tnear);
  Atn = axp_physmat_setup_mex(2,2,mu,3,1,p,np,pmodes,iside,iclosed,[1;N+1],[1;np+2],N,np+1, ...
      sx,snx,sws,swxp,tpan,Mn,[1;Mn+1],Pn,Nn,eye(3),zeros(3,1),3*Mn,nblk) ...
      + axp_physmat_setup_mex(2,4,mu,3,1,p,np,pmodes,iside,iclosed,[1;N+1],[1;np+2],N,np+1, ...
      sx,snx,sws,swxp,tpan,Mn,[1;Mn+1],Pn,Nn,eye(3),zeros(3,1),3*Mn,nblk);
  fillt = repelem(any(reshape(any(Atn~=0,2),3,Mn),1).',3,1);   % targets the master filled
  tev = Atn*sigphys; tnear(fillt) = tev(fillt);
  errT = vecnorm(reshape(tnear,3,Mn) - Tex)/max(vecnorm(Tex));
  err_trac(kk)=max(errT);
  fprintf('  near traction S''+D'' sparse corrected: max err = %.3e  (%d pts, ntcx3=%d)\n', err_trac(kk), Mn, ntcx3);
end

figure(1),clf; semilogy(np_vals,errmax,'o-k',np_vals,err_trac,'s-r');
xlabel('n_p'); ylabel('max err'); grid on; legend('near velocity','near traction','Location','best');
title(sprintf('Stokes DLPn+SLPn Neumann BVP [%s] (physmat sparse): h-refinement',shape));

function y=apply_tracfmm_close(sig, N, p, np, nang, ntcx, tcxi, ids, Ss, idd, Dd, s3dx, s3dnx, s3dw, mu, fmm_eps)
y = slpntracfmm(sig, s3dx, s3dw, [], s3dnx, fmm_eps, 1) + ...
    dlpntracfmm(sig, s3dx, s3dnx, s3dw, [], s3dnx, mu, fmm_eps, 1);
y = axps_closecorrapply_mex(3, N, p, np, nang, ntcx, tcxi, ids, Ss, sig, 1, y);
y = axps_closecorrapply_mex(3, N, p, np, nang, ntcx, tcxi, idd, Dd, sig, 1, y);
end

function U=stk(X,y,F)
d=X-y; r=vecnorm(d); U=(1/(8*pi))*(F./r+(F.'*d).*d./r.^3);
end

function T=stktrac(X,n,y,F)
d=X-y; r=vecnorm(d); rF=F.'*d; rn=sum(n.*d,1); T=-(3/(4*pi))*(rF.*rn).*d./r.^5;
end

function t=slpntracfmm(sig, s3dx, s3dw, X, NX, fmm_eps, iself)
% Naive S' traction matvec via stfmm3d.  sig is INTERLOCKED.
srcinfo.sources=s3dx; srcinfo.stoklet=reshape(sig,3,[]).*s3dw;
if iself
  U=stfmm3d(fmm_eps,srcinfo,3);       G=U.grad;     pr=U.pre;
else
  U=stfmm3d(fmm_eps,srcinfo,0,X,3);   G=U.gradtarg; pr=U.pretarg;
end
Gs=G+permute(G,[2 1 3]);
t=squeeze(sum(Gs.*reshape(NX,1,3,[]),2)) - pr.*NX;
t=t(:);
end

function t=dlpntracfmm(sig, s3dx, s3dnx, s3dw, X, NX, mu, fmm_eps, iself)
% Naive D' traction matvec via stfmm3d stresslets.  sig is INTERLOCKED.
srcinfo.sources=s3dx;
srcinfo.strslet=reshape(sig,3,[]).*s3dw;
srcinfo.strsvec=s3dnx;
if iself
  U=stfmm3d(fmm_eps,srcinfo,3);       G=U.grad;     pr=U.pre;
else
  U=stfmm3d(fmm_eps,srcinfo,0,X,3);   G=U.gradtarg; pr=U.pretarg;
end
Gs=G+permute(G,[2 1 3]);
t=mu*(squeeze(sum(Gs.*reshape(NX,1,3,[]),2)) - pr.*NX);
t=t(:);
end

function [D_ij, idxall] = closestokdlp_physmat(nt, tx, t3dx, p, np, nang, sx, snx, sws, swxp, tpan, gate, s3dx, s3dnx, s3dw, pmodes, iside, iclosed, mu, ntcx, tcxi, D_ij, idxall)
% Per-panel physical-space Stokes DLP velocity close correction.
% This mirrors axps_closestokslp_mex's pass-2 contract, but uses the existing
% offdiag physical-matrix DLP velocity path until a dedicated axps_closestokdlp_mex exists.
tx = tx(:); sx = sx(:); snx = snx(:); sws = sws(:); swxp = swxp(:); tpan = tpan(:);
N = p*np;
txr = [real(tx) imag(tx)].';
sxr = [real(sx) imag(sx)].';
kd = KDTree(txr.');
D2 = zeros(size(D_ij));
for j = 1:np
  rowsj = tcxi(j):tcxi(j+1)-1;
  jj = (j-1)*p + (1:p);
  qradii = gate*sqrt(sum(sws(jj)));
  qpoint = mean(sxr(:,jj), 2);
  idxcj = kd.ball(qpoint, qradii);
  idxcj = sort(idxcj(:));
  idxall(rowsj) = idxcj;
  T3d = t3dx(:,idxcj);
  mj = size(T3d,2);
  if mj==0, continue, end

  colj = reshape(jj(:) + (0:nang-1)*N, [], 1);
  xcsj = [zeros(2,p); imag(sx(jj)).'];
  Atj=[]; Bj=[]; Fj=[]; nrj=[];
  [Atj,Bj,Fj,nrj] = axp_offdiagphysmat_mex(2,3,3, ...
      mj, T3d, zeros(3,mj), [0;0;0], eye(3), zeros(3,1), 1, 0, ...
      nang*p, s3dx(:,colj), zeros(3,nang*p), zeros(nang*p,1), [0;0;0], eye(3), xcsj, p, ...
      sx(jj), snx(jj), sws(jj), swxp(jj), p, 1, tpan(j:j+1), 0, ...
      pmodes, iside, iclosed, mu, 1, Atj, Bj, Fj, nrj);
  Dnaive = Sto3dDLPmat_il(struct('x',T3d), struct('x',s3dx(:,colj),'nx',s3dnx(:,colj),'w',s3dw(colj)), mu);
  Dnaive(~isfinite(Dnaive)) = 0;
  D2(3*(rowsj(1)-1)+1 : 3*rowsj(end), :) = Atj - Dnaive;
end
D_ij = D2;
end
