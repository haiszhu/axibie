clearvars; format short e;
addpath('../../utils');
addpath('../../matlab');
addpath('../../external/fmm3d/matlab');   % lfmm3d/Lap3dDLPfmm

% profile clear
% profile on

% ---- geometry switch (mirror axibie/testStokesDir.m / the SLPn test) ----
% shape = 'cshape';                                   % 'sphere' | 'ellipse' | 'cshape'
shape = 'sphere';
p=16; iside=1; iclosed=0;
if strcmp(shape,'sphere')
  Z=@(t) sin(t)-1i*cos(t);                          % unit sphere
  np_vals=2:16;
  % np_vals=2:8;
  % np_vals=8:14;
  y1=[0.30;0;0.20]; q1=1.0; y2=[0.32*cos(pi/4);0.32*sin(pi/4);-0.25]; q2=-0.7;
  gv=linspace(-1.9,1.9,16);
elseif strcmp(shape,'ellipse')
  a=2; Z=@(t) a^(-1/3)*sin(t)-1i*a^(2/3)*cos(t);    % volume-conserving prolate ellipsoid
  np_vals=2:16;
  y1=[0.20;0;0.40]; q1=1.0; y2=[0.25*cos(pi/4);0.25*sin(pi/4);-0.50]; q2=-0.7;
  gv=linspace(-2.4,2.4,18);
elseif strcmp(shape,'cshape')
  lam=0.75; Z=@(t) -(1.5+cos(t)).*(-sin(lam*pi*sin(t))+1i*cos(lam*pi*sin(t)));  % c-shape (near-axis poles)
  np_vals=4:2:24;
  y1=[0.70;0;-1.20]; q1=1.0; y2=[1.30*cos(pi/4);1.30*sin(pi/4);0.20]; q2=-0.7;
  gv=linspace(-3,3,24);
end

% interior point charges -> exact exterior harmonic potential u_ex (scalar); Dirichlet data = u_ex on surface
uex=@(X) q1./(4*pi*vecnorm(X-y1)) + q2./(4*pi*vecnorm(X-y2));

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

  % 1. geometry + Dirichlet data setup
  s=[]; s.p=p; s.Z=Z; s.tpan=linspace(0,pi,np+1)'; s=quadr(s,[],'p','G');
  N=numel(s.x);
  sx=s.x(:); snx=s.nx(:); sws=s.ws(:); swxp=s.wxp(:);
  tpan=s.tpan(:); sxlo=s.Z(s.tpan(1:end-1)); sxlo=sxlo(:); sxhi=s.Z(s.tpan(2:end)); sxhi=sxhi(:);
  s3d=axisym_to_3d_quadrature(real(s.x),imag(s.x),s.ws,nang);
  phis=2*pi*(0:nang-1)/nang;                            % 3D ring UNIT normals [azimuth outer, meridian inner]
  s3dnx=[reshape(real(s.nx(:))*cos(phis),1,[]); reshape(real(s.nx(:))*sin(phis),1,[]); repmat(imag(s.nx(:)).',1,nang)];

  % % loop over patch and compute close interaction nodes, and this is essentially 2 x N computation
  % t = s; % for self evaluation, target t is source s
  % sx_r = [real(s.x(:)) imag(s.x(:))]';
  % tx_r = [real(t.x(:)) imag(t.x(:))]';
  % kd = KDTree(tx_r'); % this is something specific to axisymmetric/rotational geometry
  % % tcx_r = zeros(size(tx_r));
  % tcxi = ones(np+1,1);
  % ntcx = 0;
  % for j=1:np
  %   % targets close to this patch
  %   jj=(j-1)*p+(1:p);
  %   qradii = gate*sqrt(sum(sws(jj))); % seems to be good, maybe analysis can be done here
  %   qpoint = mean(sx_r(:,jj),2);
  %   idxcj = kd.ball( qpoint, qradii); 
  %   tcxj = tx_r(:,idxcj);
  %   %
  %   ntcxj = size(tcxj,2);
  %   tcxi(j+1) = tcxi(j) + size(tcxj,2);
  %   % tcx_r(:,tcxi_r(k):(tcxi_r(k+1)-1)) = tcxj;
  %   ntcx = ntcx + ntcxj;
  % end
  % % tcx_r = tcx_r(:,1:tcxi_r(end)-1);

  % 2. loop over patch and compute close interaction nodes, and this is essentially 2 x N computation
  % FROZEN close-correction path (uncomment to compare): t = s;                                                % self: target t is source s
  % FROZEN close-correction path (uncomment to compare): tx = t.x(:); nt = numel(tx);
  % FROZEN close-correction path (uncomment to compare): t3dx = [real(tx).'; zeros(1,nt); imag(tx).'];
  % FROZEN close-correction path (uncomment to compare): sx_r = [real(s.x(:)) imag(s.x(:))]';
  % FROZEN close-correction path (uncomment to compare): tx_r = [real(t.x(:)) imag(t.x(:))]';
  % FROZEN close-correction path (uncomment to compare): tcxi = zeros(np+1,1); ntcx = 0;                       % preallocate (INOUT, mex-style)
  % [tcxi, ntcx] = axps_closesize(nt, tx, t3dx, p, np, sx, sws, gate, tcxi, ntcx);
  % FROZEN close-correction path (uncomment to compare): [tcxi, ntcx] = axps_closesize_mex(nt, tx, t3dx, p, np, sx, sws, gate, tcxi, ntcx);

  % 3. allocate space for correction matrix, value, and I J index
  % FROZEN close-correction path (uncomment to compare): S_ij = zeros(ntcx,nang*p); % for self correction, per panel correction
  % FROZEN close-correction path (uncomment to compare): idxall = zeros(ntcx,1);

  % % actual computation
  % kd = KDTree(tx_r');
  % for j=1:np
  %   rowsj = tcxi(j):tcxi(j+1)-1;
  %   % targets close to this patch
  %   jj=(j-1)*p+(1:p);
  %   qradii = gate*sqrt(sum(sws(jj))); % seems to be good, maybe analysis can be done here
  %   qpoint = mean(sx_r(:,jj),2);
  %   idxcj = kd.ball( qpoint, qradii);
  %   idxall(tcxi(j):tcxi(j+1)-1) = idxcj;
  %   tcxj = tx_r(:,idxcj);
  %   mj = size(tcxj,2);
  %   T3d = [tcxj(1,:); zeros(1,mj); tcxj(2,:)];              % close targets, 3D at azimuth 0
  %   colj = reshape(jj(:) + (0:nang-1)*N, [], 1);           % panel j 3D source cols [azimuth outer, meridian inner]
  %   % accurate close-eval from panel j (as source) to these targets (off-diag, iforce=1)
  %   xcsj=[zeros(2,p); imag(sx(jj)).']; Atj=[];Bj=[];Fj=[];nrj=[];
  %   [Atj,Bj,Fj,nrj]=axp_offdiagphysmat_mex(1,1,1, mj, T3d, zeros(3,mj), [0;0;0],eye(3), zeros(3,1),1,0, ...
  %       nang*p, s3d.x(:,colj), zeros(3,nang*p), zeros(nang*p,1), [0;0;0],eye(3), xcsj, p, sx(jj),snx(jj),sws(jj),swxp(jj), p,1, tpan(j:j+1), 0, ...
  %       pmodes,iside,iclosed,1.0, 1, Atj,Bj,Fj,nrj);
  %   Snaive = Lap3dSLPmat(struct('x',T3d), struct('x',s3d.x(:,colj),'w',s3d.w(colj)));  % naive, same 3D sources
  %   Snaive(~isfinite(Snaive)) = 0;
  %   S_ij(rowsj,:) = Atj - Snaive;                          % accurate - naive (close correction block)
  % end

  % 4. actual computation
  % [S_ij, idxall] = axps_closeslp(nt, tx, t3dx, p, np, nang, sx, snx, sws, swxp, tpan, gate, ...
  %     s3d.x, s3dnx, s3d.w, pmodes, iside, iclosed, ntcx, tcxi, S_ij, idxall);
  % FROZEN close-correction path (uncomment to compare): [S_ij, idxall] = axps_closeslp_mex(nt, tx, t3dx, p, np, nang, sx, snx, sws, swxp, tpan, gate, ...
  % FROZEN close-correction path (uncomment to compare):     s3d.x, s3dnx, s3d.w, pmodes, iside, iclosed, ntcx, tcxi, S_ij, idxall);

  % % build actual system matrix
  % A0 = Lap3dSLPmat(s3d, s3d);
  % A0(diagind(A0)) = 0;
  % for j=1:np
  %   jj=(j-1)*p+(1:p);
  %   colj = reshape(jj(:)+(0:nang-1)*N, [], 1);
  %   rowsj = tcxi(j):tcxi(j+1)-1;
  %   idxcj  = idxall(rowsj);
  %   for c = 1:numel(idxcj)
  %     i   = idxcj(c);
  %     Pm  = reshape(S_ij(rowsj(c),:), p, nang);
  %     grows = (0:nang-1).'*N + i;
  %     Blk = zeros(nang, nang*p);
  %     for at = 1:nang
  %       Blk(at,:) = reshape(circshift(Pm,[0,at-1]), 1, []);
  %     end
  %     A0(grows, colj) = A0(grows, colj) + Blk;
  %   end
  % end
  % A = A0;
  % % keyboard

  % 5. build actual system matrix
  % FROZEN close-correction path (uncomment to compare): A = zeros(N*nang);                                    % preallocate (INOUT, mex-style)
  % A = axps_closeasm(1,1,1, nt, tx, t3dx, p, np, nang, s3d.x, s3dnx, s3d.w, 1.0, ntcx, tcxi, idxall, S_ij, 1, A);
  % FROZEN close-correction path (uncomment to compare): A = axps_closeasm_mex(1,1,1, nt, tx, t3dx, p, np, nang, s3d.x, s3dnx, s3d.w, 1.0, ntcx, tcxi, idxall, S_ij, 1, A);

  % 5'. LEVEL-2 PHYSICAL master (axp_physmat_setup_mex): assembled dense SLP self operator,
  %     LAB node-interleaved [azimuth outer, meridian inner] (== axp_physmat_mex's xA)
  nblk = N*nang;
  A = axp_physmat_setup_mex(1,1,0.0,1,1,p,np,pmodes,iside,iclosed,[1;N+1],[1;np+2],N,np+1, ...
      sx,snx,sws,swxp,tpan,0,[1;1],zeros(3,0),zeros(3,0),eye(3),zeros(3,1),nblk,nblk);

  % % 2. self operator (assembled, LAB node-interleaved) + physical-space solve  A sigma = u
  % %    (lsqminnorm: SLP meridian null space; A block-diagonalizes over modes, so == the per-mode min-norm)
  % A0 = Lap3dSLPmat(s3d, s3d);                              % naive self (accurate far, FMM-replaceable; self-diag Inf -> overwritten below)
  % A  = axp_physmat_mex(1,1,1,p,np,sx,snx,sws,swxp,tpan,sxlo,sxhi,pmodes,iside,iclosed,1.0,eye(3));   % accurate self (close-eval)
  % gate=3.75; nfill=0;
  % for j=1:np                                              % close-correct panel-by-panel: overwrite near blocks with the accurate op
  %   jj=(j-1)*p+(1:p); panlen=sum(sws(jj));                % panel k meridian nodes + arclength
  %   cl=find(abs(sx-sxlo(j))+abs(sx-sxhi(j)) < gate*panlen);           % meridian targets close to panel k (endpoint gate)
  %   trow=reshape((0:nang-1).'*N + cl(:).', [], 1);        % 3D rows: close targets, all azimuths
  %   scol=reshape((0:nang-1).'*N + jj(:).', [], 1);        % 3D cols: panel k, all azimuths
  %   A0(trow,scol)=A(trow,scol); nfill=nfill+numel(trow)*numel(scol);   % count filled entries this panel
  % end
  % A = A0;                                                 % close-corrected self operator (solve with this)
  % fprintf('  close-fill entries = %d / numel(A) = %d   (ratio %.3e)\n', nfill, numel(A), nfill/numel(A));
  
  % 6. solve
  uphys = Lap3dSLPmat(struct('x',s3d.x), struct('x',[y1 y2],'w',[1 1]))*[q1;q2];   % Dirichlet data = interior-charge potential at the 3D surface nodes ([azimuth-outer, node-inner])
  % FROZEN close-correction path (uncomment to compare): sigphys = lsqminnorm(A, uphys);
  sigphys = A\uphys;                                    % dense direct solve (backslash, as the multi dense reference)

  % % 3. 3d target eval via the physical-space off-diag operator (iforce=1 -> all targets get the block-builder eval)
  % ns=N*nang; xcs=[zeros(2,N); imag(sx).'];
  % At=[]; Bo=[]; Fo=[]; nro=[];
  % [At,Bo,Fo,nro] = axp_offdiagphysmat_mex(1,1,1, M3, P, zeros(3,M3), [0;0;0],eye(3), zeros(3,1),1,0, ...
  %     ns, s3d.x, zeros(3,ns), zeros(ns,1), [0;0;0],eye(3), xcs, N, sx,snx,sws,swxp, p,np, tpan, 0, ...
  %     pmodes,iside,iclosed,1.0, 1, At,Bo,Fo,nro);
  % u3 = (At*sigphys).';                                      % 1 x M3 exterior potential (Laplace real)

  % 7.
  % FROZEN close-correction path (uncomment to compare): ns=N*nang;
  % FROZEN close-correction path (uncomment to compare): ttx = (rho3 + 1i*z3).'; nt3 = M3;
  % FROZEN close-correction path (uncomment to compare): tcxi3 = zeros(np+1,1); ntcx3 = 0;
  % [tcxi3, ntcx3] = axps_closesize(nt3, ttx, P, p, np, sx, sws, gate, tcxi3, ntcx3);
  % FROZEN close-correction path (uncomment to compare): [tcxi3, ntcx3] = axps_closesize_mex(nt3, ttx, P, p, np, sx, sws, gate, tcxi3, ntcx3);

  % 8.
  % FROZEN close-correction path (uncomment to compare): S3 = zeros(ntcx3,nang*p); idx3 = zeros(ntcx3,1);
  % [S3, idx3] = axps_closeslp(nt3, ttx, P, p, np, nang, sx, snx, sws, swxp, tpan, gate, s3d.x, s3dnx, s3d.w, pmodes, iside, iclosed, ntcx3, tcxi3, S3, idx3);
  % FROZEN close-correction path (uncomment to compare): [S3, idx3] = axps_closeslp_mex(nt3, ttx, P, p, np, nang, sx, snx, sws, swxp, tpan, gate, s3d.x, s3dnx, s3d.w, pmodes, iside, iclosed, ntcx3, tcxi3, S3, idx3);

  % % 9.
  % At = zeros(nt3, ns);
  % At = axps_closeasm(1,1,1, nt3, ttx, P, p, np, nang, s3d.x, s3dnx, s3d.w, 1.0, ntcx3, tcxi3, idx3, S3, 0, At);
  % u3 = (At*sigphys).';

  % 9'. matrix-free
  % FROZEN close-correction path (uncomment to compare): u3v = zeros(nt3,1);
  % u3v = axps_closeopasm(1,1,1, nt3, ttx, P, p, np, nang, s3d.x, s3dnx, s3d.w, ntcx3, tcxi3, idx3, S3, sigphys, 0, 1e-14, u3v);
  % FROZEN close-correction path (uncomment to compare): u3v = Lap3dSLPfmm(struct('x',P), struct('x',s3d.x,'w',s3d.w), sigphys, 1e-14);
  % FROZEN close-correction path (uncomment to compare): u3v = axps_closecorrapply_mex(1, nt3, p, np, nang, ntcx3, tcxi3, idx3, S3, sigphys, 0, u3v);
  % FROZEN close-correction path (uncomment to compare): u3 = u3v.';

  % 9''. LEVEL-2 PHYSICAL master eval (iinter=3): master fills its internal near zone (far rows
  %      stay 0) -> dense naive base at all grid targets, master rows replace it (multi-reference
  %      eval pattern: naive full + near close replacement)
  Aev = axp_physmat_setup_mex(1,1,0.0,3,1,p,np,pmodes,iside,iclosed,[1;N+1],[1;np+2],N,np+1, ...
      sx,snx,sws,swxp,tpan,M3,[1;M3+1],P,zeros(3,M3),eye(3),zeros(3,1),M3,nblk);
  u3v = Lap3dSLPmat(struct('x',P), struct('x',s3d.x,'w',s3d.w))*sigphys;   % dense naive base (accurate where master is far)
  fillm = any(Aev~=0,2);                                % targets the master filled (its internal near zone)
  uev = Aev*sigphys; u3v(fillm) = uev(fillm);
  u3 = u3v.';

  % 5. log10 error
  Uex3=uex(P); err3=abs(u3-Uex3)/max(abs(Uex3));
  inb=d3<0.1; errmax(kk)=max(err3);
  fprintf('Laplace SLP all-modes Dirichlet BVP [%s] (delta mex): N=%d, np=%d, pmodes=%d, M3=%d\n',shape,N,np,pmodes,M3);
  if any(inb), fprintf('  near band (d<0.1):  max err = %.3e  (%d pts)\n',max(err3(inb)),nnz(inb)); end
  fprintf('  far field (d>=0.1): max err = %.3e  (%d pts)\n',max(err3(~inb)),nnz(~inb));
end

figure(1),clf; semilogy(np_vals,errmax,'o-k'); xlabel('n_p'); ylabel('max err'); grid on;
title(sprintf('Laplace SLP all-modes Dirichlet BVP [%s]: h-refinement, p_{modes}=2 n_p',shape));

% scatter3 log10 error (last refinement)
figure(2),clf; hold on;
thm=linspace(0,2*pi,49);
Lx=real(s.x)*cos(thm); Ly=real(s.x)*sin(thm); Lz=imag(s.x)*ones(size(thm));
mesh(Lx,Ly,Lz,'EdgeColor',[0.55 0.55 0.55],'FaceColor','none','EdgeAlpha',0.35);
scatter3(P(1,:),P(2,:),P(3,:),18,log10(max(err3,1e-17)),'filled');
plot3(y1(1),y1(2),y1(3),'p','MarkerSize',16,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
plot3(y2(1),y2(2),y2(3),'p','MarkerSize',16,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
cb=colorbar; cb.Label.String='log_{10}|u_h-u_{exact}|/max|u_{exact}|';
clim([-16 -12]); axis equal; view(35,18); grid on; colormap('jet');
xlabel('x'); ylabel('y'); zlabel('z'); title(sprintf('Laplace SLP all-modes [%s]: 3D target-grid potential error',shape));

% profile viewer

% cmp = getPyPlot_cMap('rainbow', [], [], '/Users/hzhu/.venvs/mpl/bin/python'); 
% colormap(cmp)
% cmp = getPyPlot_cMap('rainbow', [], [], '"/Users/hzhu/.pyenv/versions/3.11.13/bin/python"');
% colormap(cmp)

% exportgraphics(figure(1),'axissymslap_lap_slp_convergence2.png','Resolution',200)
% exportgraphics(figure(2),'axissymslap_lap_slp_error2.png','Resolution',200)
