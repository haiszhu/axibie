clearvars; format short e;
addpath('../../../axibie/utils');
addpath('../../utils');
addpath('../../matlab');
addpath('../../external/kdtree/toolbox')
addpath('../../external/fmm3d/matlab');

% ---- geometry switch (mirror the stok_dlp per-panel test) ----
shape = 'cshape';                                   % 'sphere' | 'ellipse' | 'cshape'
shape = 'sphere';
p=16; iside=1; iclosed=0;
if strcmp(shape,'sphere')
  Z=@(t) sin(t)-1i*cos(t);                          % unit sphere
  np_vals=2:16;
  y1=[0.30;0;0.20]; y2=[0.32*cos(pi/4);0.32*sin(pi/4);-0.25];
  gv=linspace(-1.9,1.9,16);
elseif strcmp(shape,'ellipse')
  a=2; Z=@(t) a^(-1/3)*sin(t)-1i*a^(2/3)*cos(t);    % volume-conserving prolate ellipsoid
  np_vals=2:16;
  y1=[0.20;0;0.40]; y2=[0.25*cos(pi/4);0.25*sin(pi/4);-0.50];
  gv=linspace(-2.4,2.4,18);
elseif strcmp(shape,'cshape')
  lam=0.75; Z=@(t) -(1.5+cos(t)).*(-sin(lam*pi*sin(t))+1i*cos(lam*pi*sin(t)));  % c-shape (near-axis poles)
  np_vals=4:2:24;
  y1=[0.70;0;-1.20]; y2=[1.30*cos(pi/4);1.30*sin(pi/4);0.20];
  gv=linspace(-3,3,24);
end
F1=(y2-y1)/norm(y2-y1); F2=-F1;                     % force- & torque-free Stokeslet pair: F2=-F1, F1 || (y2-y1)

% interior Stokeslets -> exact exterior Stokes velocity u_ex; Neumann data = traction on surface.
% single-layer ansatz u = S[sigma]; solve the traction Neumann system (-1/2 I + S') sigma = t.
% (-1/2 I + S') = adjoint of the interior-Dirichlet DLP (-1/2 I + D) -> 6-dim rigid-motion
% deficiency: solvable only for data with ZERO net force AND torque (hence the F2=-F1
% collinear pair above); sigma is non-unique, lsqminnorm picks a particular solution
% (any particular sigma gives the same exterior field u = S[sigma]).
uex=@(X) stk(X,y1,F1)+stk(X,y2,F2);

% ---- exterior 3D target grid (once): outside the revolved meridian, >1e-2 off the surface ----
mer=Z(linspace(0,pi,400).');
[GX,GY,GZ]=meshgrid(gv,gv,gv); rho3=hypot(GX(:),GY(:)); z3=GZ(:); zt=rho3+1i*z3;
dmer=inf(size(zt)); for j=1:numel(mer), dmer=min(dmer,abs(zt-mer(j))); end
ext3=(~inpolygon(rho3,z3,real(mer),imag(mer))) & (dmer>1e-2);
P=[GX(ext3).';GY(ext3).';GZ(ext3).']; M3=size(P,2);
rho3=hypot(P(1,:),P(2,:)); z3=P(3,:); phi3=atan2(P(2,:),P(1,:)); d3=dmer(ext3).';

errmax=nan(1,numel(np_vals));           % near velocity (post-eval of u = S[sigma])
err_trac=nan(1,numel(np_vals));         % traction referee (S'[sigma] at random 3D normals)
gate = 2.0;
for kk=1:numel(np_vals)
  np=np_vals(kk); pmodes=2*np; nmodes=2*pmodes+1; nang=nmodes; mu=1;

  % 1. geometry + Neumann (LAB traction) data setup
  s=[]; s.p=p; s.Z=Z; s.tpan=linspace(0,pi,np+1)'; s=quadr(s,[],'p','G');
  N=numel(s.x);
  sx=s.x(:); snx=s.nx(:); sws=s.ws(:); swxp=s.wxp(:);
  tpan=s.tpan(:); sxlo=s.Z(s.tpan(1:end-1)); sxlo=sxlo(:); sxhi=s.Z(s.tpan(2:end)); sxhi=sxhi(:);
  s3d=axisym_to_3d_quadrature(real(s.x),imag(s.x),s.ws,nang);
  phis=2*pi*(0:nang-1)/nang;                            % 3D ring UNIT normals [azimuth outer, meridian inner]
  s3dnx=[reshape(real(s.nx(:))*cos(phis),1,[]); reshape(real(s.nx(:))*sin(phis),1,[]); repmat(imag(s.nx(:)).',1,nang)];
  nrm=real(snx).'; nzm=imag(snx).';
  tlab=zeros(3,N*nang);                                 % exact surface traction (LAB, INTERLOCKED node-inner azimuth-outer)
  for i=1:nang
    th=2*pi*(i-1)/nang; X=s3d.x(:,(i-1)*N+(1:N)); N3=[nrm*cos(th); nrm*sin(th); nzm];
    tlab(:,(i-1)*N+(1:N)) = stktrac(X,N3,y1,F1)+stktrac(X,N3,y2,F2);
  end
  tphys=tlab(:);

  % 2. per-panel SELF system matrix (-1/2 I + S'): axps_naivestoksdlpn_physmat
  %    (theta=0 block-rows with own-panel target rows excluded BY INDEX, block-circulant
  %    tile with the rotz(phi_a) 3x3 conjugation inside), then per-panel close-overwrite
  %    (axps_closestoksdlpn_panel_mex, exterior one-sided limit -> the -1/2 I jump rides
  %    in S'): wrapper evaluated ONCE at theta=0, written into ALL nang azimuth row-blocks
  %    with the circulant column shift AND the rotz conjugation W_a = R(phi_a) W0 R(phi_a)^T
  % LEVEL-2 master reference (uncomment to compare):
  % nblk = 3*N*nang;
  % A = axp_physmat_setup_mex(2,2,mu,1,1,p,np,pmodes,iside,iclosed,[1;N+1],[1;np+2],N,np+1, ...
  %     sx,snx,sws,swxp,tpan,0,[1;1],zeros(3,0),zeros(3,0),eye(3),zeros(3,1),nblk,nblk);   % -1/2 I + S'
  phib = 2*pi*(0:nang-1)/nang;
  t0 = [real(sx).'; zeros(1,N); imag(sx).'];            % theta=0 meridian 3D coords
  n0 = [real(snx).'; zeros(1,N); imag(snx).'];          % theta=0 3D surface normals
  [Aselfs, Aselfd] = axps_naivestoksdlpn_physmat_mex(N, sx, t0, n0, p, np, nang, sx, snx, sws, pmodes, mu, 1, 3*N*nang, [], []);
  A = Aselfs;                                           % S'-only: ignore Aselfd
  for j = 1:np
    idxj = (j-1)*p + (1:p);
    Tj = find(abs(sx - sxlo(j)) + abs(sx - sxhi(j)) < 1.85*sum(sws(idxj)));
    if isempty(Tj), continue, end
    icl = 0; if j==1, icl = -1; elseif j==np, icl = 1; end
    P0 = [real(sx(Tj)).'; zeros(1,numel(Tj)); imag(sx(Tj)).'];   % theta=0 3D coords
    N0 = [real(snx(Tj)).'; zeros(1,numel(Tj)); imag(snx(Tj)).']; % theta=0 3D surface normals
    skx = sx(idxj);
    [Sn3,Dn3] = axps_closestoksdlpn_panel_mex(numel(Tj), sx(Tj).', P0, N0, p, nang, skx, ...
                snx(idxj), sws(idxj), swxp(idxj), sxlo(j), sxhi(j), [],[],[],[],[], pmodes, iside, icl, mu, ...
                [],[],[],1,[],[]);   % iform=1: ON-surface self rows -> pole graded sweep OFF (blockmat mesh)
    W0 = Sn3;                                           % S'-only: ignore Dn3
    rowT = reshape((3*(Tj(:)-1) + (1:3)).', [], 1);     % target xyz rows (theta=0 block)
    cbj  = reshape((3*(idxj(:)-1) + (1:3)).', [], 1);   % panel xyz cols inside an azimuth block
    for a = 1:nang                                      % evaluated once, written nang times
      ca = cos(phib(a)); sa = sin(phib(a));             % rotz conjugation W_a = R W0 R^T
      Wx = W0(1:3:end,:); Wy = W0(2:3:end,:);           % target rows: rotate xy by +phi_a
      Wa = W0; Wa(1:3:end,:) = ca*Wx - sa*Wy; Wa(2:3:end,:) = sa*Wx + ca*Wy;
      Vx = Wa(:,1:3:end); Vy = Wa(:,2:3:end);           % source cols: right-multiply by R^T
      Wa(:,1:3:end) = ca*Vx - sa*Vy; Wa(:,2:3:end) = sa*Vx + ca*Vy;
      bglob = mod((a-1)+(0:nang-1), nang);              % circulant azimuth-block shift
      cols = reshape(3*N*bglob + cbj, [], 1);
      A((a-1)*3*N + rowT, cols) = Wa;
    end
  end

  % 3. solve the traction Neumann system (lsqminnorm: (-1/2 I + S') has the 6-dim
  %    rigid-motion nullspace; any particular sigma gives the same field u = S[sigma])
  sigphys = lsqminnorm(A, tphys);

  % 4. per-panel post-eval of the VELOCITY u = S[sigma] (explicit naive modal matrix +
  %    per-panel close-overwrite via axps_closestoksdlp_panel, routed by panel position)
  ttx = (rho3 + 1i*z3).';
  [Astargs, Adtargs] = axps_naivestoksdlp_physmat_mex(M3, ttx, P, p, np, nang, sx, snx, sws, pmodes, mu, 0, 3*M3, [], []);
  Atargs = Astargs;                                    % SLP-only: ignore Adtargs
  for j = 1:np
    idxj = (j-1)*p + (1:p);
    Tj = find(abs(ttx - sxlo(j)) + abs(ttx - sxhi(j)) < 1.85*sum(sws(idxj)));
    if isempty(Tj), continue, end
    icl = 0; if j==1, icl = -1; elseif j==np, icl = 1; end
    skx = sx(idxj);
    [S3p,D3p] = axps_closestoksdlp_panel_mex(numel(Tj), ttx(Tj).', P(:,Tj), [], p, nang, skx, ...
                [],[],[], sxlo(j), sxhi(j), [],[],[],[],[], pmodes, iside, icl, mu, ...
                [],[],[],[],[],[]);
    ring = reshape(idxj(:) + (0:nang-1)*N, [], 1);
    rowi = reshape((3*(Tj(:)-1) + (1:3)).', [], 1);
    coli = reshape((3*(ring-1)  + (1:3)).', [], 1);
    Atargs(rowi, coli) = S3p;
  end
  u3p = reshape(Atargs*sigphys, 3, M3);

  % 5. traction referee: t = S'[sigma] at RANDOM 3D target normals (n_theta ~= 0),
  %    explicit naive traction modal matrix + per-panel close-overwrite via
  %    axps_closestoksdlpn_panel.  Exact reference = interior-Stokeslet surface traction.
  rng(11); Nr = randn(3,M3); Nr = Nr./vecnorm(Nr);
  [Asn, Adn] = axps_naivestoksdlpn_physmat_mex(M3, ttx, P, Nr, p, np, nang, sx, snx, sws, pmodes, mu, 0, 3*M3, [], []);
  Atn = Asn;                                           % S'-only: ignore Adn
  for j = 1:np
    idxj = (j-1)*p + (1:p);
    Tj = find(abs(ttx - sxlo(j)) + abs(ttx - sxhi(j)) < 1.85*sum(sws(idxj)));
    if isempty(Tj), continue, end
    icl = 0; if j==1, icl = -1; elseif j==np, icl = 1; end
    skx = sx(idxj);
    [Sn3,Dn3] = axps_closestoksdlpn_panel_mex(numel(Tj), ttx(Tj).', P(:,Tj), Nr(:,Tj), p, nang, skx, ...
                snx(idxj), sws(idxj), swxp(idxj), sxlo(j), sxhi(j), [],[],[],[],[], pmodes, iside, icl, mu, ...
                [],[],[],0,[],[]);
    ring = reshape(idxj(:) + (0:nang-1)*N, [], 1);
    rowi = reshape((3*(Tj(:)-1) + (1:3)).', [], 1);
    coli = reshape((3*(ring-1)  + (1:3)).', [], 1);
    Atn(rowi, coli) = Sn3;
  end
  t3p = reshape(Atn*sigphys, 3, M3);
  Tex = stktrac(P, Nr, y1, F1) + stktrac(P, Nr, y2, F2);

  % 6. errors
  Uex3 = uex(P);
  fprintf('Stokes (-1/2 I + S'') traction Neumann BVP [%s] (per-panel): N=%d, np=%d, pmodes=%d, M3=%d\n',shape,N,np,pmodes,M3);
  err3p = vecnorm(u3p-Uex3)/max(vecnorm(Uex3));
  errmax(kk)=max(err3p);
  fprintf('  velocity  (S) eval over ALL %d targets: max err = %.3e\n', M3, max(err3p));
  errT = vecnorm(t3p-Tex)/max(vecnorm(Tex));
  err_trac(kk)=max(errT);
  fprintf('  traction (S'') referee, random 3D n over ALL %d targets: max err = %.3e\n', M3, max(errT));

  % per-loop scatter3 log10 error (which region is causing issue)
  figure(3); clf,
  subplot(1,2,1)
  thm=linspace(0,2*pi,49);
  mesh(real(s.x)*cos(thm), real(s.x)*sin(thm), imag(s.x)*ones(size(thm)), ...
       'EdgeColor',[0.55 0.55 0.55],'FaceColor','none','EdgeAlpha',0.35);
  scatter3(P(1,:),P(2,:),P(3,:),8,log10(max(err3p,1e-17)),'filled');
  axis equal; view(35,18); grid on; colormap('jet'); colorbar;
  caxis([-16 0])
  title(sprintf('S'' Neumann np=%d (max %.1e)', np, max(err3p)));
  drawnow;

  % per-loop scatter3 log10 traction error (random 3D normals)
  figure(3);
  subplot(1,2,2)
  mesh(real(s.x)*cos(thm), real(s.x)*sin(thm), imag(s.x)*ones(size(thm)), ...
       'EdgeColor',[0.55 0.55 0.55],'FaceColor','none','EdgeAlpha',0.35); hold on
  scatter3(P(1,:),P(2,:),P(3,:),8,log10(max(errT,1e-17)),'filled');
  axis equal; view(35,18); grid on; colormap('jet'); colorbar;
  caxis([-16 0])
  title(sprintf('traction (random n) np=%d (max %.1e)', np, max(errT)));
  drawnow;
end

figure(1),clf; semilogy(np_vals,errmax,'o-k',np_vals,err_trac,'s-r');
xlabel('n_p'); ylabel('max err'); grid on; legend('velocity (S)','traction (S'')','Location','best');
title(sprintf('Stokes (-1/2 I + S'') traction Neumann BVP [%s] (per-panel): h-refinement',shape));

function U=stk(X,y,F)
d=X-y; r=vecnorm(d); U=(1/(8*pi))*(F./r+(F.'*d).*d./r.^3);
end

function T=stktrac(X,n,y,F)
d=X-y; r=vecnorm(d); rF=F.'*d; rn=sum(n.*d,1); T=-(3/(4*pi))*(rF.*rn).*d./r.^5;
end
