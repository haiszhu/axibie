clearvars; format short e;
addpath('../../../axibie/utils');
addpath('../../utils');
addpath('../../matlab');
addpath('../../external/kdtree/toolbox')
addpath('../../external/fmm3d/matlab');

% ---- geometry switch (mirror axibie/testStokesDir.m / the SLP & SLPn tests) ----
shape = 'cshape';                                   % 'sphere' | 'ellipse' | 'cshape'
shape = 'sphere';
p=16; iside=1; iclosed=0;
if strcmp(shape,'sphere')
  Z=@(t) sin(t)-1i*cos(t);                          % unit sphere
  np_vals=2:16;
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

% interior point charges -> exact exterior harmonic potential u_ex; Neumann data g = dn u_ex.
% Combined-field Neumann: represent u = (S + D)[sigma], solve (-1/2 I + S' + D') sigma = dn u_ex.
% Per-panel physical build: S' carries the exterior jump (-1/2 I) via the one-sided limit;
% D' is continuous; the hypersingular nullspace D'[1] = 0 is enforced by a row-sum
% diagonal correction on the assembled physical D' self matrix (m=0 blockmat fix analogue).
uex  =@(X) q1./(4*pi*vecnorm(X-y1)) + q2./(4*pi*vecnorm(X-y2));
gradu=@(X) -q1*(X-y1)./(4*pi*vecnorm(X-y1).^3) - q2*(X-y2)./(4*pi*vecnorm(X-y2).^3);

% ---- exterior 3D target grid (once): outside the revolved meridian, >1e-2 off the surface ----
mer=Z(linspace(0,pi,400).');                        % fine meridian for inpolygon + distance
[GX,GY,GZ]=meshgrid(gv,gv,gv); rho3=hypot(GX(:),GY(:)); z3=GZ(:); zt=rho3+1i*z3;
dmer=inf(size(zt)); for j=1:numel(mer), dmer=min(dmer,abs(zt-mer(j))); end
ext3=(~inpolygon(rho3,z3,real(mer),imag(mer))) & (dmer>1e-2);
P=[GX(ext3).';GY(ext3).';GZ(ext3).']; M3=size(P,2);
rho3=hypot(P(1,:),P(2,:)); z3=P(3,:); phi3=atan2(P(2,:),P(1,:)); d3=dmer(ext3).';

errmax=nan(1,numel(np_vals));
for kk=1:numel(np_vals)
  np=np_vals(kk); pmodes=2*np; nmodes=2*pmodes+1; nang=nmodes;

  % 1. geometry + Neumann (dn u_ex) data at the 3D surface nodes
  s=[]; s.p=p; s.Z=Z; s.tpan=linspace(0,pi,np+1)'; s=quadr(s,[],'p','G');
  N=numel(s.x);
  sx=s.x(:); snx=s.nx(:); sws=s.ws(:); swxp=s.wxp(:);
  tpan=s.tpan(:); sxlo=s.xlo(:); sxhi=s.xhi(:);
  s3d=axisym_to_3d_quadrature(real(s.x),imag(s.x),s.ws,nang);
  phis=2*pi*(0:nang-1)/nang;                            % 3D ring UNIT normals [azimuth outer, meridian inner]
  s3dnx=[reshape(real(s.nx(:))*cos(phis),1,[]); reshape(real(s.nx(:))*sin(phis),1,[]); repmat(imag(s.nx(:)).',1,nang)];
  gphys = sum(s3dnx.*gradu(s3d.x),1).';                 % Neumann data ([azimuth-outer, node-inner])

  % 2. per-panel SELF system matrix (-1/2 I + S' + D'): naive S'/D' physical matrices
  %    in ONE call (theta=0 block-row + circulant tile inside), then per-panel
  %    close-overwrite (axps_closelapsdlpn_panel) written into ALL azimuth row-blocks
  %    with the circulant shift.  S' and D' kept separate for the D' nullspace fix.
  t0 = [real(sx).'; zeros(1,N); imag(sx).'];            % theta=0 meridian 3D coords
  n0 = [real(snx).'; zeros(1,N); imag(snx).'];          % theta=0 3D surface normals
  [Aselfs, Aselfd] = axps_naivelapsdlpn_physmat_mex(N, sx, t0, n0, p, np, nang, sx, snx, sws, pmodes, 1, N*nang, [], []);
  for j = 1:np
    idxj = (j-1)*p + (1:p);
    Tj = find(abs(sx - sxlo(j)) + abs(sx - sxhi(j)) < 1.85*sum(sws(idxj)));
    if isempty(Tj), continue, end
    icl = 0; if j==1, icl = -1; elseif j==np, icl = 1; end
    P0 = [real(sx(Tj)).'; zeros(1,numel(Tj)); imag(sx(Tj)).'];   % theta=0 3D coords
    N0 = [real(snx(Tj)).'; zeros(1,numel(Tj)); imag(snx(Tj)).']; % theta=0 3D surface normals
    skx = sx(idxj);
    [Sn3,Dn3] = axps_closelapsdlpn_panel_mex(numel(Tj), sx(Tj).', P0, N0, p, nang, skx, ...
                [],[],[], sxlo(j), sxhi(j), [],[],[],[],[], pmodes, iside, icl, ...
                [],[],[],[],[],[]);
    for a = 1:nang                                      % evaluated once, written nang times
      rows = (a-1)*N + Tj;
      cols = reshape(idxj(:) + mod((a-1)+(0:nang-1), nang)*N, 1, []);
      Aselfs(rows, cols) = Sn3;
      Aselfd(rows, cols) = Dn3;
    end
  end
  % hypersingular nullspace: D'[1] = 0 -> subtract each row-sum from the diagonal
  rs = sum(Aselfd, 2);
  Aselfd(1:N*nang+1:end) = Aselfd(1:N*nang+1:end) - rs.';
  Aself = Aselfs + Aselfd;                    % = -1/2 I + S' + D' (S' exterior limit carries the jump)

  % 3. solve
  sigphys = Aself\gphys;

  % 4. per-panel post-eval of the combined potential u = (S + D)[sigma]
  ttx = (rho3 + 1i*z3).';
  [Astargs, Adtargs] = axps_naivelapsdlp_physmat_mex(M3, ttx, P, p, np, nang, sx, snx, sws, pmodes, 0, M3, [], []);
  Atargs = Astargs + Adtargs;
  for j = 1:np
    idxj = (j-1)*p + (1:p);
    Tj = find(abs(ttx - sxlo(j)) + abs(ttx - sxhi(j)) < 1.85*sum(sws(idxj)));
    if isempty(Tj), continue, end
    icl = 0; if j==1, icl = -1; elseif j==np, icl = 1; end
    skx = sx(idxj);
    [S3p,D3p] = axps_closelapsdlp_panel_mex(numel(Tj), ttx(Tj).', P(:,Tj), [], p, nang, skx, ...
                [],[],[], sxlo(j), sxhi(j), [],[],[],[],[], pmodes, iside, icl, ...
                [],[],[],[],[],[]);
    ring = reshape(idxj(:) + (0:nang-1)*N, [], 1);
    Atargs(Tj, ring) = S3p + D3p;
  end
  u3p = (Atargs*sigphys).';

  % 4'. traction referee at field targets with RANDOM 3D normals (n_theta ~= 0):
  %     n.grad u = (S'+D')[sigma], exact reference from the point charges.
  %     First live exercise of the azimuthal-normal cross-term in BOTH the naive
  %     physmat (ifself=0) and the close wrapper (t3dnx = random normals).
  rng(11); Nr = randn(3,M3); Nr = Nr./vecnorm(Nr);      % one random unit normal per target
  [Asn_t, Adn_t] = axps_naivelapsdlpn_physmat_mex(M3, ttx, P, Nr, p, np, nang, sx, snx, sws, pmodes, 0, M3, [], []);
  Atn = Asn_t + Adn_t;
  for j = 1:np
    idxj = (j-1)*p + (1:p);
    Tj = find(abs(ttx - sxlo(j)) + abs(ttx - sxhi(j)) < 1.85*sum(sws(idxj)));
    if isempty(Tj), continue, end
    icl = 0; if j==1, icl = -1; elseif j==np, icl = 1; end
    skx = sx(idxj);
    [Sn3,Dn3] = axps_closelapsdlpn_panel_mex(numel(Tj), ttx(Tj).', P(:,Tj), Nr(:,Tj), p, nang, skx, ...
                [],[],[], sxlo(j), sxhi(j), [],[],[],[],[], pmodes, iside, icl, ...
                [],[],[],[],[],[]);
    ring = reshape(idxj(:) + (0:nang-1)*N, [], 1);
    Atn(Tj, ring) = Sn3 + Dn3;
  end
  gn3p = (Atn*sigphys).';
  gnex = sum(Nr.*gradu(P),1);                           % exact n.grad u_ex (analytic)

  % 5. log10 error
  Uex3=uex(P);
  fprintf('Laplace combined (S''+D'') all-modes Neumann BVP [%s] (per-panel): N=%d, np=%d, pmodes=%d, M3=%d\n',shape,N,np,pmodes,M3);
  err3p=abs(u3p-Uex3)/max(abs(Uex3));
  errmax(kk)=max(err3p);
  fprintf('  per-panel eval (all-axi): max err over ALL %d targets = %.3e\n', M3, max(errmax(kk)));
  errn3=abs(gn3p-gnex)/max(abs(gnex));
  fprintf('  4'' traction, random 3D normals: max err over ALL %d targets = %.3e\n', M3, max(errn3));

  % per-loop scatter3 log10 error (which region is causing issue)
  figure(3); clf,
  subplot(1,2,1)
  thm=linspace(0,2*pi,49);
  mesh(real(s.x)*cos(thm), real(s.x)*sin(thm), imag(s.x)*ones(size(thm)), ...
       'EdgeColor',[0.55 0.55 0.55],'FaceColor','none','EdgeAlpha',0.35);
  scatter3(P(1,:),P(2,:),P(3,:),8,log10(max(err3p,1e-17)),'filled');
  axis equal; view(35,18); grid on; colormap('jet'); colorbar;
  caxis([-16 0])
  title(sprintf('S''+D'' Neumann np=%d (max %.1e)', np, max(err3p)));
  drawnow;

  % per-loop scatter3 log10 traction error (random 3D normals)
  figure(3); 
  subplot(1,2,2)
  mesh(real(s.x)*cos(thm), real(s.x)*sin(thm), imag(s.x)*ones(size(thm)), ...
       'EdgeColor',[0.55 0.55 0.55],'FaceColor','none','EdgeAlpha',0.35); hold on
  scatter3(P(1,:),P(2,:),P(3,:),8,log10(max(errn3,1e-17)),'filled');
  axis equal; view(35,18); grid on; colormap('jet'); colorbar;
  caxis([-16 0])
  title(sprintf('traction (random n) np=%d (max %.1e)', np, max(errn3)));
  drawnow;
end

figure(1),clf; semilogy(np_vals,errmax,'o-k'); xlabel('n_p'); ylabel('max err'); grid on;
title(sprintf('Laplace combined (S''+D'') all-modes Neumann BVP [%s]: h-refinement, p_{modes}=2 n_p',shape));

% scatter3 log10 error (last refinement)
figure(2),clf; hold on;
thm=linspace(0,2*pi,49);
Lx=real(s.x)*cos(thm); Ly=real(s.x)*sin(thm); Lz=imag(s.x)*ones(size(thm));
mesh(Lx,Ly,Lz,'EdgeColor',[0.55 0.55 0.55],'FaceColor','none','EdgeAlpha',0.35);
scatter3(P(1,:),P(2,:),P(3,:),18,log10(max(err3p,1e-17)),'filled');
plot3(y1(1),y1(2),y1(3),'p','MarkerSize',16,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
plot3(y2(1),y2(2),y2(3),'p','MarkerSize',16,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
cb=colorbar; cb.Label.String='log_{10}|u_h-u_{exact}|/max|u_{exact}|';
clim([-16 -12]); axis equal; view(35,18); grid on; colormap('jet');
xlabel('x'); ylabel('y'); zlabel('z'); title(sprintf('Laplace combined (S''+D'') all-modes Neumann [%s]: 3D target-grid potential error',shape));
ylim([-2 3])

% exportgraphics(figure(1),'axissymslap_lap_dlpn_convergence_panel.png','Resolution',200)
% exportgraphics(figure(2),'axissymslap_lap_dlpn_error_panel.png','Resolution',200)
