clearvars; format short e;
addpath('../../utils');
addpath('../../matlab');

% ---- geometry switch (mirror test_axissymslap_lap_dlp_bvp_physmat.m) ----
shape = 'sphere';                                   % 'sphere' | 'ellipse' | 'cshape'
p=16; iside=1; iclosed=0; mu=1;
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
  lam=0.75; Z=@(t) -(1.5+cos(t)).*(-sin(lam*pi*sin(t))+1i*cos(lam*pi*sin(t)));  % c-shape (near-axis poles)
  np_vals=6:2:24;
  y1=[0.70;0;-1.20]; y2=[1.30*cos(pi/4);1.30*sin(pi/4);0.20];
  gv=linspace(-3,3,24);
end
F1=[1;-0.7;0.5]; F2=[0.4;0.9;-0.6];                 % interior Stokeslet forces (shape-independent)

% interior Stokeslets -> exact exterior Stokes velocity u_ex; Neumann data = traction on surface.
uex=@(X) stk(X,y1,F1)+stk(X,y2,F2);

% ---- exterior 3D target grid (once): outside the revolved meridian, >1e-2 off the surface ----
mer=Z(linspace(0,pi,400).');                        % fine meridian for inpolygon + distance
[GX,GY,GZ]=meshgrid(gv,gv,gv); rho3=hypot(GX(:),GY(:)); z3=GZ(:); zt=rho3+1i*z3;
dmer=inf(size(zt)); for j=1:numel(mer), dmer=min(dmer,abs(zt-mer(j))); end
ext3=(~inpolygon(rho3,z3,real(mer),imag(mer))) & (dmer>1e-2);
P=[GX(ext3).';GY(ext3).';GZ(ext3).']; M3=size(P,2);
rho3=hypot(P(1,:),P(2,:)); z3=P(3,:); phi3=atan2(P(2,:),P(1,:)); d3=dmer(ext3).';

errmax=nan(1,numel(np_vals));
gate = 2.0;
for kk=1:numel(np_vals)
  np=np_vals(kk); pmodes=2*np; nmodes=2*pmodes+1; nang=nmodes;

  % 1. geometry + Neumann data setup (LAB traction at ALL ring nodes, INTERLOCKED)
  s=[]; s.p=p; s.Z=Z; s.tpan=linspace(0,pi,np+1)'; s=quadr(s,[],'p','G');
  N=numel(s.x);
  sx=s.x(:); snx=s.nx(:); sws=s.ws(:); swxp=s.wxp(:); tpan=s.tpan(:);
  s3d=axisym_to_3d_quadrature(real(s.x),imag(s.x),s.ws,nang); nr=real(s.nx).'; nz=imag(s.nx).';
  phis=2*pi*(0:nang-1)/nang;                            % 3D ring UNIT normals [azimuth outer, meridian inner]
  s3dnx=[reshape(real(snx)*cos(phis),1,[]); reshape(real(snx)*sin(phis),1,[]); repmat(imag(snx).',1,nang)];
  nsf=N*nang;
  tlab=zeros(3,nsf);
  for i=1:nang
    th=2*pi*(i-1)/nang; X=s3d.x(:,(i-1)*N+(1:N)); N3=[nr*cos(th);nr*sin(th);nz];
    tlab(:,(i-1)*N+(1:N)) = stktrac(X,N3,y1,F1)+stktrac(X,N3,y2,F2);
  end
  tphys=tlab(:);

  % 2. loop over patch and compute close interaction nodes (self: azimuth-0 orbit reps)
  % FROZEN close-correction path (uncomment to compare): tx = sx; nt = N;
  % FROZEN close-correction path (uncomment to compare): t3dx  = [real(sx).'; zeros(1,N); imag(sx).'];
  % FROZEN close-correction path (uncomment to compare): tn3dx = [real(snx).'; zeros(1,N); imag(snx).'];      % target normals = surface normals (n_theta=0 exact)
  % FROZEN close-correction path (uncomment to compare): tcxi = zeros(np+1,1); ntcx = 0;                       % preallocate (INOUT, mex-style)
  % FROZEN close-correction path (uncomment to compare): [tcxi, ntcx] = axps_closesize_mex(nt, tx, t3dx, p, np, sx, sws, gate, tcxi, ntcx);

  % 3. allocate space for correction matrix, value, and I J index
  % FROZEN close-correction path (uncomment to compare): S_ij = zeros(3*ntcx, 3*nang*p);                       % nc=3 INTERLOCKED close-correction blocks
  % FROZEN close-correction path (uncomment to compare): idxall = zeros(ntcx,1);

  % 4. actual computation (S' traction close corrections)
  % [S_ij, idxall] = axps_closestokslpn(nt, tx, t3dx, tn3dx, p, np, nang, sx, snx, sws, swxp, tpan, gate, ...
  %     s3d.x, s3dnx, s3d.w, pmodes, iside, iclosed, mu, ntcx, tcxi, S_ij, idxall);
  % FROZEN close-correction path (uncomment to compare): [S_ij, idxall] = axps_closestokslpn_mex(nt, tx, t3dx, tn3dx, p, np, nang, sx, snx, sws, swxp, tpan, gate, ...
  % FROZEN close-correction path (uncomment to compare):     s3d.x, s3dnx, s3d.w, pmodes, iside, iclosed, mu, ntcx, tcxi, S_ij, idxall);

  % 5. matrix-free operator: stfmm3d naive S' at the sources (self term excluded by the FMM)
  %    + circulant corrapply (nc=3: azimuth shift + R(phi_a) block conjugation).  The operator
  %    IS the exterior limit -1/2 I + S' (close-eval AT the nodes, iside=1 -- jump inside).
  % FROZEN close-correction path (uncomment to compare): applyA = @(x) axps_closecorrapply_mex(3, N, p, np, nang, ntcx, tcxi, idxall, S_ij, x, 1, ...
  % FROZEN close-correction path (uncomment to compare):                   slpntracfmm(x, s3d.x, s3d.w, [], s3dnx, 1e-14, 1));

  % 5'. LEVEL-2 PHYSICAL master (axp_physmat_setup_mex): assembled dense SELF operator, LAB
  %     node-interleaved (== axp_physmat_mex's xA); the exterior jump rides inside (-1/2 I + S')
  nblk = 3*N*nang;
  A = axp_physmat_setup_mex(2,2,mu,1,1,p,np,pmodes,iside,iclosed,[1;N+1],[1;np+2],N,np+1, ...
      sx,snx,sws,swxp,tpan,0,[1;1],zeros(3,0),zeros(3,0),eye(3),zeros(3,1),nblk,nblk);

  % 6. solve (dense min-norm, as the multi dense reference: pure S' is rank-deficient, the
  %    surface-normal density is the SLP pressure-gauge null vector)
  % FROZEN close-correction path (uncomment to compare): nh = s3dnx(:)/norm(s3dnx(:));
  % FROZEN close-correction path (uncomment to compare): [sigphys,flag,relres,iter] = gmres(@(x) applyA(x) + nh*(nh'*x), tphys, [], 1e-12, 200);
  % FROZEN close-correction path (uncomment to compare): fprintf('  gmres: flag=%d  iters=%d  relres=%.3e  sig.n=%.2e\n', flag, iter(2), relres, nh'*sigphys);
  sigphys = lsqminnorm(A, tphys);

  % 7. kdtree near-target detection ONCE (ball per panel) -- the SAME near set feeds the
  %    velocity AND traction verifications; far grid targets get pure FMM, no correction.
  %    (build the grid tree AFTER the axps_* self passes: those clobber the global kdtree)
  ttx = (rho3 + 1i*z3).';
  axk_kdtree_build_mex(1, [rho3(:), z3(:), zeros(M3,1)], 64);
  isnear = false(M3,1);
  for k = 1:np
    jj = (k-1)*p + (1:p);
    qr = mean(real(s.x(jj))); qz = mean(imag(s.x(jj)));
    qradii = gate*sqrt(sum(s.ws(jj)));                  % == the closesize detection rule: naive
    idxcin = axk_kdtree_ball_mex([qr, qz, 0], qradii, M3);  % needs Nv~34/acosh(chi) -> d ~ sqrt(h)
    isnear(idxcin) = true;
  end
  Pn = P(:,isnear); Mn = nnz(isnear); ttn = ttx(isnear);
  rng(0); Nn = randn(3,Mn); Nn = Nn./vecnorm(Nn);       % random target normals (traction check)

  % 8. close sizes ONCE at the near set + BOTH compute passes (velocity S3, traction Stn)
  % FROZEN close-correction path (uncomment to compare): tcxi3 = zeros(np+1,1); ntcx3 = 0;
  % FROZEN close-correction path (uncomment to compare): [tcxi3, ntcx3] = axps_closesize_mex(Mn, ttn, Pn, p, np, sx, sws, gate, tcxi3, ntcx3);
  % FROZEN close-correction path (uncomment to compare): S3 = zeros(3*ntcx3, 3*nang*p); idx3 = zeros(ntcx3,1);
  % FROZEN close-correction path (uncomment to compare): [S3, idx3] = axps_closestokslp_mex(Mn, ttn, Pn, p, np, nang, sx, snx, sws, swxp, tpan, gate, ...
  % FROZEN close-correction path (uncomment to compare):     s3d.x, s3dnx, s3d.w, pmodes, iside, iclosed, mu, ntcx3, tcxi3, S3, idx3);
  % FROZEN close-correction path (uncomment to compare): Stn = zeros(3*ntcx3, 3*nang*p); itn = zeros(ntcx3,1);
  % FROZEN close-correction path (uncomment to compare): [Stn, itn] = axps_closestokslpn_mex(Mn, ttn, Pn, Nn, p, np, nang, sx, snx, sws, swxp, tpan, gate, ...
  % FROZEN close-correction path (uncomment to compare):     s3d.x, s3dnx, s3d.w, pmodes, iside, iclosed, mu, ntcx3, tcxi3, Stn, itn);

  % 9'. matrix-free velocity eval: Sto3dSLPfmm at ALL targets + LEVEL-2 master eval (iinter=3) at
  %     the kdtree near set: the master fills its internal near zone (far rows stay 0), those rows
  %     replace the FMM values (naive full + near replacement, the multi-reference eval pattern)
  sigblk = reshape(reshape(sigphys,3,[]).',[],1);       % interlocked -> component-block
  u3blk  = Sto3dSLPfmm(struct('x',P), struct('x',s3d.x,'w',s3d.w), sigblk, 1e-14);
  u3v    = reshape(reshape(u3blk,[],3).',[],1);         % component-block -> interlocked (3*M3)
  rown = reshape(3*(find(isnear).'-1)+(1:3)',[],1);     % near-set rows of the global u
  un = u3v(rown);
  % FROZEN close-correction path (uncomment to compare): un = axps_closecorrapply_mex(3, Mn, p, np, nang, ntcx3, tcxi3, idx3, S3, sigphys, 0, un);
  Aev = axp_physmat_setup_mex(2,1,mu,3,1,p,np,pmodes,iside,iclosed,[1;N+1],[1;np+2],N,np+1, ...
      sx,snx,sws,swxp,tpan,Mn,[1;Mn+1],Pn,zeros(3,Mn),eye(3),zeros(3,1),3*Mn,nblk);
  fillm = repelem(any(reshape(any(Aev~=0,2),3,Mn),1).',3,1);   % targets the master filled (interlocked rows)
  uev = Aev*sigphys; un(fillm) = uev(fillm);
  u3v(rown) = un;
  U3 = reshape(u3v, 3, M3);

  % 10. log10 velocity error (near = kdtree set, corrected; far = pure FMM)
  Uex3=uex(P); err3=vecnorm(U3-Uex3)/max(vecnorm(Uex3));
  errmax(kk)=max(err3(isnear));
  fprintf('Stokes SLPn all-modes Neumann BVP [%s] (physmat, matrix-free): N=%d, np=%d, pmodes=%d, M3=%d\n',shape,N,np,pmodes,M3);
  fprintf('  kdtree-near (corrected): max err = %.3e  (%d pts)\n', max(err3(isnear)), Mn);
  fprintf('  far (pure FMM):          max err = %.3e  (%d pts)\n', max(err3(~isnear)), nnz(~isnear));

  % 11. S' traction verification at the SAME near set (naive FMM base + master eval replacement)
  Tex = stktrac(Pn,Nn,y1,F1) + stktrac(Pn,Nn,y2,F2);
  t3v = slpntracfmm(sigphys, s3d.x, s3d.w, Pn, Nn, 1e-14, 0);
  % FROZEN close-correction path (uncomment to compare): t3v = axps_closecorrapply_mex(3, Mn, p, np, nang, ntcx3, tcxi3, itn, Stn, sigphys, 0, t3v);
  Atn = axp_physmat_setup_mex(2,2,mu,3,1,p,np,pmodes,iside,iclosed,[1;N+1],[1;np+2],N,np+1, ...
      sx,snx,sws,swxp,tpan,Mn,[1;Mn+1],Pn,Nn,eye(3),zeros(3,1),3*Mn,nblk);
  fillt = repelem(any(reshape(any(Atn~=0,2),3,Mn),1).',3,1);   % targets the master filled
  tev = Atn*sigphys; t3v(fillt) = tev(fillt);
  errT = vecnorm(reshape(t3v,3,Mn) - Tex)/max(vecnorm(Tex));
  fprintf('  S'' traction at near set: max err = %.3e  (%d pts)\n', max(errT), Mn);
end

figure(1),clf; semilogy(np_vals,errmax,'o-k'); xlabel('n_p'); ylabel('max err'); grid on;
title(sprintf('Stokes SLPn all-modes Neumann BVP [%s] (physmat): h-refinement, p_{modes}=2 n_p',shape));

% scatter3 log10 error (last refinement)
figure(2),clf; hold on;
thm=linspace(0,2*pi,49);
Lx=real(s.x)*cos(thm); Ly=real(s.x)*sin(thm); Lz=imag(s.x)*ones(size(thm));
mesh(Lx,Ly,Lz,'EdgeColor',[0.55 0.55 0.55],'FaceColor','none','EdgeAlpha',0.35);
scatter3(P(1,:),P(2,:),P(3,:),18,log10(max(err3,1e-17)),'filled');
plot3(y1(1),y1(2),y1(3),'p','MarkerSize',16,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
plot3(y2(1),y2(2),y2(3),'p','MarkerSize',16,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
cb=colorbar; cb.Label.String='log_{10}|u_h-u_{exact}|/max|u_{exact}|';
clim([-16 -8]); axis equal; view(35,18); grid on; colormap('jet');
xlabel('x'); ylabel('y'); zlabel('z'); title(sprintf('Stokes SLPn all-modes [%s]: 3D target-grid velocity error',shape));

% exportgraphics(figure(1),'axissymsstok_stok_slpn_convergence2.png','Resolution',200)
% exportgraphics(figure(2),'axissymsstok_stok_slpn_error2.png','Resolution',200)

function U=stk(X,y,F)
d=X-y; r=vecnorm(d); U=(1/(8*pi))*(F./r+(F.'*d).*d./r.^3);
end
function T=stktrac(X,n,y,F)
d=X-y; r=vecnorm(d); rF=F.'*d; rn=sum(n.*d,1); T=-(3/(4*pi))*(rF.*rn).*d./r.^5;
end
function t=slpntracfmm(sig, s3dx, s3dw, X, NX, fmm_eps, iself)
% naive S' traction matvec via stfmm3d:  t = -p n + (grad u + grad u^T) n   (stfmm3d mu=1;
% S' is mu-independent).  iself=1: at the SOURCES themselves (self term excluded by the FMM,
% == the naive operator's zeroed diagonal); iself=0: at targets X.  sig INTERLOCKED (3*ns).
srcinfo.sources=s3dx; srcinfo.stoklet=reshape(sig,3,[]).*s3dw;
if iself
  U=stfmm3d(fmm_eps,srcinfo,3);       G=U.grad;     pr=U.pre;
else
  U=stfmm3d(fmm_eps,srcinfo,0,X,3);   G=U.gradtarg; pr=U.pretarg;
end
Gs=G+permute(G,[2 1 3]);                              % only the symmetric part enters
t=squeeze(sum(Gs.*reshape(NX,1,3,[]),2)) - pr.*NX;
t=t(:);
end
