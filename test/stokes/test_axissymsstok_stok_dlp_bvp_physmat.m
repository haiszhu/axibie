clearvars; format short e;
addpath('../../../axibie/utils');
addpath('../../utils');
addpath('../../matlab');
addpath('../../external/kdtree/toolbox')
addpath('../../external/fmm3d/matlab');   % lfmm3d/Lap3dDLPfmm

% ---- geometry switch (mirror axibie/testStokesDir.m / the SLP & SLPn tests) ----
shape = 'cshape';                                   % 'sphere' | 'ellipse' | 'cshape'
shape = 'sphere';
p=16; iside=1; iclosed=0;
if strcmp(shape,'sphere')
  Z=@(t) sin(t)-1i*cos(t);                          % unit sphere
  np_vals=2:14;
  y1=[0.30;0;0.20]; q1=1.0; y2=[0.32*cos(pi/4);0.32*sin(pi/4);-0.25]; q2=-0.7;
  gv=linspace(-1.9,1.9,16);
elseif strcmp(shape,'ellipse')
  a=2; Z=@(t) a^(-1/3)*sin(t)-1i*a^(2/3)*cos(t);    % volume-conserving prolate ellipsoid
  np_vals=2:14;
  y1=[0.20;0;0.40]; q1=1.0; y2=[0.25*cos(pi/4);0.25*sin(pi/4);-0.50]; q2=-0.7;
  gv=linspace(-2.4,2.4,18);
elseif strcmp(shape,'cshape')
  lam=0.75; Z=@(t) -(1.5+cos(t)).*(-sin(lam*pi*sin(t))+1i*cos(lam*pi*sin(t)));  % c-shape (near-axis poles)
  np_vals=4:2:24;
  y1=[0.70;0;-1.20]; q1=1.0; y2=[1.30*cos(pi/4);1.30*sin(pi/4);0.20]; q2=-0.7;
  gv=linspace(-3,3,24);
end

% interior Stokeslets -> exact exterior Stokes velocity u_ex; Dirichlet data = u_ex on surface.
% (combined field u = (D_+ + S)[sigma], full rank; mirror test_axissymsstok_stok_dlp_bvp.m)
F1=[1;-0.7;0.5]; F2=[0.4;0.9;-0.6];                 % interior Stokeslet forces (mu=1)
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
  np=np_vals(kk); pmodes=2*np; nmodes=2*pmodes+1; nang=nmodes; mu=1;

  % 1. geometry + Dirichlet data setup
  s=[]; s.p=p; s.Z=Z; s.tpan=linspace(0,pi,np+1)'; s=quadr(s,[],'p','G');
  N=numel(s.x);
  sx=s.x(:); snx=s.nx(:); sws=s.ws(:); swxp=s.wxp(:);
  tpan=s.tpan(:); sxlo=s.Z(s.tpan(1:end-1)); sxlo=sxlo(:); sxhi=s.Z(s.tpan(2:end)); sxhi=sxhi(:);
  s3d=axisym_to_3d_quadrature(real(s.x),imag(s.x),s.ws,nang);
  phis=2*pi*(0:nang-1)/nang;                            % 3D ring UNIT normals [azimuth outer, meridian inner]
  s3dnx=[reshape(real(s.nx(:))*cos(phis),1,[]); reshape(real(s.nx(:))*sin(phis),1,[]); repmat(imag(s.nx(:)).',1,nang)];

  % 2. per-panel SELF system matrix (D_+ + S): axps_naivestoksdlp_physmat (theta=0
  %    block-rows with own-panel target rows excluded BY INDEX, block-circulant tile
  %    with the rotz(phi_a) 3x3 conjugation inside), then per-panel close-overwrite
  %    (axps_closestoksdlp_panel_mex, exterior one-sided limit -> the +1/2 I DLP jump
  %    rides inside): wrapper evaluated ONCE at theta=0, written into ALL nang azimuth
  %    row-blocks with the circulant column shift AND the rotz conjugation
  %    W_a = R(phi_a) W0 R(phi_a)^T (VECTOR kernel; scalar Laplace write had no twist)
  % LEVEL-2 master reference (uncomment to compare):
  % nblk = 3*N*nang;
  % As = axp_physmat_setup_mex(2,1,mu,1,1,p,np,pmodes,iside,iclosed,[1;N+1],[1;np+2],N,np+1, ...
  %     sx,snx,sws,swxp,tpan,0,[1;1],zeros(3,0),zeros(3,0),eye(3),zeros(3,1),nblk,nblk);   % S
  % Ad = axp_physmat_setup_mex(2,3,mu,1,1,p,np,pmodes,iside,iclosed,[1;N+1],[1;np+2],N,np+1, ...
  %     sx,snx,sws,swxp,tpan,0,[1;1],zeros(3,0),zeros(3,0),eye(3),zeros(3,1),nblk,nblk);   % +1/2 I + D
  % Aself = As + Ad;
  phib = 2*pi*(0:nang-1)/nang;
  t0 = [real(sx).'; zeros(1,N); imag(sx).'];            % theta=0 meridian 3D coords
  [Aselfs, Aselfd] = axps_naivestoksdlp_physmat_mex(N, sx, t0, p, np, nang, sx, snx, sws, pmodes, mu, 1, 3*N*nang, [], []);
  Aself = Aselfs + Aselfd;
  for j = 1:np
    idxj = (j-1)*p + (1:p);
    Tj = find(abs(sx - sxlo(j)) + abs(sx - sxhi(j)) < 1.85*sum(sws(idxj)));
    if isempty(Tj), continue, end
    iclosed = 0; if j==1, iclosed = -1; elseif j==np, iclosed = 1; end
    P0 = [real(sx(Tj)).'; zeros(1,numel(Tj)); imag(sx(Tj)).'];   % theta=0 3D coords
    skx = sx(idxj);
    [S3p,D3p] = axps_closestoksdlp_panel_mex(numel(Tj), sx(Tj).', P0, [], p, nang, skx, ...
                [],[],[], sxlo(j), sxhi(j), [],[],[],[],[], pmodes, iside, iclosed, mu, ...
                [],[],[],1,[],[]);   % iform=1: ON-surface self rows -> pole graded sweep OFF (blockmat mesh)
    W0 = S3p + D3p;
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
      Aself((a-1)*3*N + rowT, cols) = Wa;
    end
  end

  % 6. solve (dense direct backslash, as the multi dense reference: combined (D_+ + S) full-rank)
  Ulab = uex(s3d.x);                       % Dirichlet data = Stokeslet velocity at the 3D surface nodes
  uphys = Ulab(:);                         % LAB xyz innermost, node inner, azimuth outer (== self rows)
  sigphys = Aself\uphys;

  % 3. per-panel post-eval (ALL-axi quadrature): ONE full modal blockmat at every
  %       target, then per panel OVERWRITE the blockmat's own panel block with the
  %       per-panel wrapper on that panel's 1.85 near zone (exact swap: the subtracted
  %       block is a slice of the same full-blockmat matrices); wrapper routed by
  %       panel position (j=1 -> iclosed=-1, j=np -> +1, else 0)
  % axp_physmat_setup_mex & axp_physmat_setup_mex for 9''
  %%%%% naive system matrix in physical space (axps_naivestoksdlp_physmat, ifself=0:
  %%%%% whole-particle kernel call per mode, fold + rotations accumulated per (md,b))
  ttx = (rho3 + 1i*z3).';
  [Ats, Atd] = axps_naivestoksdlp_physmat_mex(M3, ttx, P, p, np, nang, sx, snx, sws, pmodes, mu, 0, 3*M3, [], []);
  Atargs = Ats + Atd;
  %%%%%
  for j = 1:np
    idxj = (j-1)*p + (1:p);
    Tj = find(abs(ttx - sxlo(j)) + abs(ttx - sxhi(j)) < 1.85*sum(sws(idxj)));
    if isempty(Tj), continue, end
    iclosed = 0; if j==1, iclosed = -1; elseif j==np, iclosed = 1; end
    skx = sx(idxj);
    [S3p,D3p] = axps_closestoksdlp_panel_mex(numel(Tj), ttx(Tj).', P(:,Tj), [], p, nang, skx, ...
                [],[],[], sxlo(j), sxhi(j), [],[],[],[],[], pmodes, iside, iclosed, mu, ...
                [],[],[],[],[],[]);
    ring = reshape(idxj(:) + (0:nang-1)*N, [], 1);      % ring nodes (node inner, azimuth outer)
    rowi = reshape((3*(Tj(:)-1) + (1:3)).', [], 1);     % target xyz rows (component innermost)
    coli = reshape((3*(ring-1)  + (1:3)).', [], 1);     % ring xyz cols (== wrapper col order)
    Atargs(rowi, coli) = S3p + D3p;
  end
  u3p = reshape(Atargs*sigphys, 3, M3);

  % 4. log10 error
  Uex3=uex(P);
  fprintf('Stokes combined (D_+ + S) all-modes Dirichlet BVP [%s] (per-panel): N=%d, np=%d, pmodes=%d, M3=%d\n',shape,N,np,pmodes,M3);
  err3p=vecnorm(u3p-Uex3)/max(vecnorm(Uex3));
  errmax(kk)=max(err3p);
  fprintf('  4  per-panel eval (all-axi): max err over ALL %d targets = %.3e\n', M3, max(errmax(kk)));

  % per-loop scatter3 log10 error of the 9''' solution (which region is causing issue)
  figure(3); clf,
  thm=linspace(0,2*pi,49);
  mesh(real(s.x)*cos(thm), real(s.x)*sin(thm), imag(s.x)*ones(size(thm)), ...
       'EdgeColor',[0.55 0.55 0.55],'FaceColor','none','EdgeAlpha',0.35);
  scatter3(P(1,:),P(2,:),P(3,:),8,log10(max(err3p,1e-17)),'filled');
  axis equal; view(35,18); grid on; colormap('jet'); colorbar;
  title(sprintf('9'''''' np=%d (max %.1e)', np, max(err3p)));
  drawnow;
end

figure(1),clf; semilogy(np_vals,errmax,'o-k'); xlabel('n_p'); ylabel('max err'); grid on;
title(sprintf('Stokes combined (D_+ + S) all-modes Dirichlet BVP [%s]: h-refinement, p_{modes}=2 n_p',shape));

% scatter3 log10 error (last refinement)
figure(2),clf; hold on;
thm=linspace(0,2*pi,49);
Lx=real(s.x)*cos(thm); Ly=real(s.x)*sin(thm); Lz=imag(s.x)*ones(size(thm));
mesh(Lx,Ly,Lz,'EdgeColor',[0.55 0.55 0.55],'FaceColor','none','EdgeAlpha',0.35);
scatter3(P(1,:),P(2,:),P(3,:),18,log10(max(err3p,1e-17)),'filled');
plot3(y1(1),y1(2),y1(3),'p','MarkerSize',16,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
plot3(y2(1),y2(2),y2(3),'p','MarkerSize',16,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
cb=colorbar; cb.Label.String='log_{10}|u_h-u_{exact}|/max|u_{exact}|';
clim([-15 -8]); axis equal; view(35,18); grid on; colormap('jet');
xlabel('x'); ylabel('y'); zlabel('z'); title(sprintf('Stokes combined (D_+ + S) all-modes [%s]: 3D target-grid velocity error',shape));
ylim([-2 3])

% cmp = getPyPlot_cMap('rainbow', [], [], '/Users/hzhu/.venvs/mpl/bin/python'); 
% colormap(cmp)
% cmp = getPyPlot_cMap('rainbow', [], [], '"/Users/hzhu/.pyenv/versions/3.11.13/bin/python"');
% colormap(cmp)

% exportgraphics(figure(1),'axissymstok_stok_dlp_convergence.png','Resolution',200)
% exportgraphics(figure(2),'axissymstok_stok_dlp_error.png','Resolution',200)

function U=stk(X,y,F)
d=X-y; r=vecnorm(d); U=(1/(8*pi))*(F./r+(F.'*d).*d./r.^3);
end