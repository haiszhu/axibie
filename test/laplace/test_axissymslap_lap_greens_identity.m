addpath('../../matlab/')
addpath('../../utils/')
addpath('../../external/fmm3d/matlab/')

clear all

% set parameters
p    = 16;
% np   = 8; % 7.54e-09
% np   = 12; % 5.12e-11
np   = 16; % 1.02e-13
M    = 4*np;
nphi = 2*M+1;
type = 'torus';
%type = 'sphere';
side = 'e';
N    = np*p;
k    = 5; isclosed = 0;
iside = double(side=='e');
%k    = 1; isclosed = -1;
if strcmp(type,'sphere'), vel = -45; else, vel = 55; end   % sphere targets sit under the south cap

% generate torus/sphere geometry (full 3d + axisymmetric)
so.a = 1.0; so.b = 0.5;
s = [];
switch type
  case 'torus'
    s.p = p; s.Z  = @(t) (so.a+so.b*cos(t-pi/2)) + 1i*(so.b*sin(t-pi/2)); s.tpan = linspace(0, 2*pi, np+1)'; 

    % define density function
    [f, gradf]= torus_surf_den(so,side); % for GRF
    phi = inoutfun_torus(so,side);
  case 'sphere'
    s.p = p; s.Z  = @(t) sin(t) - 1i*cos(t); s.tpan = linspace(0, pi, np+1)';

    % define density function
    [f, gradf]= sphere_surf_den(side); % for GRF
    phi = inoutfun_sphere(side);
end
s   = quadr(s,[],'p','G');
s3d = axisym_to_3d_quadrature( real(s.x), imag(s.x), s.ws, real(s.nx), imag(s.nx), nphi);

% all on surface node
sx = s3d.x; snx = s3d.nx; sw = s3d.w;

sigma = sum(snx.*gradf(sx),1)'; tau = -f(sx)';     % densities for GRF (col vecs)

% exterior GRF test
nt = 50000;
rng(2);
switch type
  case 'torus'
    c0 = [0.8520;0.4919;0.4819];
  case 'sphere'
    tm = (k-0.5)*pi/np; c0 = [sin(tm); 0; -cos(tm)];
end
t3d.x = c0 + 2.0*(rand(3,nt)-1/2); % random targets, shrink 2*sqrt(sum(pan{k}.w)), or increase nt to observe closer targets...
inout_flag = phi(t3d.x(1,:),t3d.x(2,:),t3d.x(3,:));
t3d_inout.x = t3d.x(:,inout_flag);
X = reshape(s3d.x(1,:), numel(s.x), nphi); Y = reshape(s3d.x(2,:), numel(s.x), nphi); Z = reshape(s3d.x(3,:), numel(s.x), nphi); X=X([1:end 1],[1:end 1]); Y=Y([1:end 1],[1:end 1]); Z=Z([1:end 1],[1:end 1]);

% solution
fhom = f(t3d_inout.x)';

% naive evaluation
eps = 1e-14;
slpval = Lap3dSLPfmm(t3d_inout,s3d,sigma,eps); 
dlpval = Lap3dDLPfmm(t3d_inout,s3d,tau,eps);
u3d_inout = slpval + dlpval; 

% reorder FFT output to match
modes = -M:M;
sigma_hat = fft(reshape(sigma,N,nphi)',2*M+1,1)/nphi;
sigma_hat = fftshift(sigma_hat, 1);  % puts negative freqs first
tau_hat   = fft(reshape(tau,N,nphi)',2*M+1,1)/nphi;
tau_hat   = fftshift(tau_hat, 1);

% error visualization
if side == 'i', err3d = abs(u3d_inout - fhom); end
if side == 'e', err3d = abs(u3d_inout + fhom); end
fprintf('N=%d (num of source), interior GRF test at %d close target pt: u err = %.3g\n', size(s3d.x,2), numel(t3d_inout.x)/3, max(err3d));
figure(1),clf,
subplot(1,3,1)
scatter3(t3d_inout.x(1,:),t3d_inout.x(2,:),t3d_inout.x(3,:),7,log10(err3d),'filled'),colorbar,
hold on, h=surf(X,Y,Z); set(h,'ambientstrength',0.7,'facelighting','gouraud'); shading interp; set(h,'FaceColor','w'); lightangle(45,0); axis equal vis3d; set(gca,'clipping','off'); xlabel('x'); ylabel('y'); zlabel('z'); colorbar;
caxis([-16 2]),title(sprintf('naive full-3D (max err = %.1e)', max(err3d)))
view(150,vel)   % close-target arc on near side, unoccluded (sphere: from below the south cap)
colorbar off; ax=gca; cb=colorbar; drawnow; ax.Position=ax.Position; p=cb.Position; cb.Position=[p(1) p(2)+p(4)/4 p(3) p(4)/2];
cmp = getPyPlot_cMap('rainbow', [], [], '"/Users/hzhu/.pyenv/versions/3.11.13/bin/python"'); colormap(cmp)

% process target, s3d = 3d_to_axisym_targets
t_inout.x     = hypot(t3d_inout.x(1,:), t3d_inout.x(2,:)).' + 1i*t3d_inout.x(3,:).';
t_inout.theta = atan2(t3d_inout.x(2,:), t3d_inout.x(1,:)).'; % angle

% classify targets into points close to k-th panel 
idx = (k-1)*s.p + (1:s.p);
ik  = abs(t_inout.x - s.xlo(k)) + abs(t_inout.x - s.xhi(k)) < 1.85*sum(s.ws(idx));
sqn_idx = find(ik);  no_sqn_idx = find(~ik);

%
tk_inout.x     = t_inout.x(sqn_idx);
tk_inout.theta = t_inout.theta(sqn_idx);
tk3d_inout.x   = t3d_inout.x(:,sqn_idx);

% skip contribution from k-th panel
nk = numel(tk_inout.x);
uk_skip_k = zeros(nk,1);
As = axls_slp_blockmat_nmode_mex(nk, tk_inout.x(:).', s.p, s.np, s.x, s.nx, s.ws, s.wxp, s.tpan, s.xlo, s.xhi, M, iside, 1, [], 1, k);
Ad = axls_dlp_blockmat_nmode_mex(nk, tk_inout.x(:).', s.p, s.np, s.x, s.nx, s.ws, s.wxp, s.tpan, s.xlo, s.xhi, M, iside, 1, [], 1, k);
for jj = 1:numel(modes)
  n = modes(jj);
  uk_skip_k = uk_skip_k + (As(:,:,abs(n)+1)*sigma_hat(jj,:).' + Ad(:,:,abs(n)+1)*tau_hat(jj,:).') .* exp(1i*n*tk_inout.theta);
end
uk_skip_k = real(uk_skip_k);

% add back naive k-th panel contribution
uk_only_k = zeros(nk,1);
idxk = (k-1)*s.p + (1:s.p);
sk = []; sk.x = s.x(idxk); sk.nx = s.nx(idxk); sk.ws = s.ws(idxk);
for jj = 1:numel(modes)
  n = modes(jj);
  Sn_k = LapSLPAxiMattmp(tk_inout, sk, n);
  Dn_k = LapDLPAxiMattmp(tk_inout, sk, n);
  utmp_k = Sn_k*sigma_hat(jj,idxk).' + Dn_k*tau_hat(jj,idxk).';
  uk_only_k = uk_only_k + utmp_k.*exp(1i*n*tk_inout.theta);
end
uk_only_k = real(uk_only_k);
u2 = uk_skip_k + uk_only_k;
err2 = abs(u2 + fhom(sqn_idx));
fprintf('without k-th panel close-eval corrected: max err at %d close targets = %.3g\n', nk, max(err2(:)));
figure(1),
subplot(1,3,2)
scatter3(tk3d_inout.x(1,:), tk3d_inout.x(2,:), tk3d_inout.x(3,:), 12, log10(err2), 'filled'); colorbar, colormap(cmp)
hold on, h=surf(X,Y,Z); set(h,'ambientstrength',0.7,'facelighting','gouraud'); shading interp; set(h,'FaceColor','w'); lightangle(45,0); axis equal vis3d; set(gca,'clipping','off'); xlabel('x'); ylabel('y'); zlabel('z'); colorbar;
axis equal vis3d; set(gca,'clipping','off'); xlabel('x'); ylabel('y'); zlabel('z'); caxis([-16 2]);
title(sprintf('naive eval (panel k=%d / %d)', k, np))
view(150,vel)
colorbar off; ax=gca; cb=colorbar; drawnow; ax.Position=ax.Position; p=cb.Position; cb.Position=[p(1) p(2)+p(4)/4 p(3) p(4)/2];

% use close k-th panel contribution instead
skx = s.x((k-1)*s.p + (1:s.p));
[As3d, Ad3d] = axps_closelapsdlp_panel(nk, tk_inout.x, tk3d_inout.x, [], s.p, nphi, skx, [], [], [], ...
                                          s.xlo(k), s.xhi(k), [], [], [], [], [], M, iside, isclosed, ...
                                          [], [], [], [], [], []);
% [As3d, Ad3d] = axps_closelapsdlp_panel(nk, tk_inout.x, tk3d_inout.x, [], s.p, nphi, skx, [], [], [], ...
%                                           s.xlo(k), s.xhi(k), [], [], [], [], [], M, iside, isclosed, ...
%                                           [], [], [], [], [], []);
cols = (k-1)*s.p;
ring = reshape((cols+(1:s.p)).' + (0:nphi-1)*N, [], 1);
uk_only_k_close = As3d*sigma(ring) + Ad3d*tau(ring);
u3 = uk_skip_k + uk_only_k_close;
err3 = abs(u3 + fhom(sqn_idx));
fprintf('with k-th panel close-eval corrected: max err at %d close targets = %.3g\n', nk, max(err3(:)));
figure(1),
subplot(1,3,3)
scatter3(tk3d_inout.x(1,:), tk3d_inout.x(2,:), tk3d_inout.x(3,:), 12, log10(err3), 'filled'); colorbar, colormap(cmp)
hold on, h=surf(X,Y,Z); set(h,'ambientstrength',0.7,'facelighting','gouraud'); shading interp; set(h,'FaceColor','w'); lightangle(45,0); axis equal vis3d; set(gca,'clipping','off'); xlabel('x'); ylabel('y'); zlabel('z'); colorbar;
axis equal vis3d; set(gca,'clipping','off'); xlabel('x'); ylabel('y'); zlabel('z'); caxis([-16 2]);
title(sprintf('with close eval (panel k=%d / %d)', k, np))
view(150,vel)
colorbar off; ax=gca; cb=colorbar; drawnow; ax.Position=ax.Position; p=cb.Position; cb.Position=[p(1) p(2)+p(4)/4 p(3) p(4)/2];

% keyboard
% exportgraphics(figure(1),'axissymslap_lap_greens_identity.png','Resolution',200)
% exportgraphics(figure(1),'axissymslap_lap_greens_identity_pole.png','Resolution',200)

p = s.p;                                 % plotting code overwrote p with cb.Position
nk = numel(sqn_idx);
Pk = tk3d_inout.x;
nblk = N*nphi;
AevS = axp_physmat_setup_mex(1,1,0.0,3,1,p,np,M,iside,0,[1;N+1],[1;np+2],N,np+1, ...
    s.x(:),s.nx(:),s.ws(:),s.wxp(:),s.tpan(:),nk,[1;nk+1],Pk,zeros(3,nk),eye(3),zeros(3,1),nk,nblk);
AevD = axp_physmat_setup_mex(1,3,0.0,3,1,p,np,M,iside,0,[1;N+1],[1;np+2],N,np+1, ...
    s.x(:),s.nx(:),s.ws(:),s.wxp(:),s.tpan(:),nk,[1;nk+1],Pk,zeros(3,nk),eye(3),zeros(3,1),nk,nblk);
u4 = AevS*sigma + AevD*tau;
fillm = any(AevS~=0,2) | any(AevD~=0,2);
u4(~fillm) = u3d_inout(sqn_idx(~fillm));
err4 = abs(u4 + fhom(sqn_idx));
fprintf('level-2 master at %d close targets: max err = %.3g (%d rows filled)\n', nk, max(err4), nnz(fillm));

nk = numel(sqn_idx(170));
Pk = tk3d_inout.x(:,170);
nblk = N*nphi;
AevS = axp_physmat_setup_mex(1,1,0.0,3,1,p,np,M,iside,0,[1;N+1],[1;np+2],N,np+1, ...
    s.x(:),s.nx(:),s.ws(:),s.wxp(:),s.tpan(:),nk,[1;nk+1],Pk,zeros(3,nk),eye(3),zeros(3,1),nk,nblk);
AevD = axp_physmat_setup_mex(1,3,0.0,3,1,p,np,M,iside,0,[1;N+1],[1;np+2],N,np+1, ...
    s.x(:),s.nx(:),s.ws(:),s.wxp(:),s.tpan(:),nk,[1;nk+1],Pk,zeros(3,nk),eye(3),zeros(3,1),nk,nblk);
u5 = AevS*sigma + AevD*tau;
err5 = abs(u5 + fhom(sqn_idx(170)))

function G = LapSLPAxiMattmp(t,s,n)
n = abs(n);                                        % kernel is even in the mode index
x1 = real(t.x(:)); x2 = imag(t.x(:));              % target (column)
y1 = real(s.x(:))'; y2 = imag(s.x(:))';            % source (row)
r2  = bsxfun(@minus,x1,y1).^2 + bsxfun(@minus,x2,y2).^2;
chi = 1 + r2./(2*bsxfun(@times,x1,y1));
tmp = 2./(chi+1); [K,E] = ellipke(tmp); sm = sqrt(tmp);
Qa = sm.*K;                                      % Q_{-1/2}
Qb = chi.*Qa - (chi+1).*sm.*E;                   % Q_{1/2}
for k = 1:n
  Qc = (2*k*chi.*Qb - (k-0.5).*Qa)./(k+0.5); Qa = Qb; Qb = Qc;
end
VK = Qa;
S  = 2*sqrt(bsxfun(@rdivide,y1,x1));
G  = 1/(4*pi)*(VK.*S)*diag(s.ws);
end

function A = LapDLPAxiMattmp(t,s,n)
n = abs(n);                                        % kernel is even in the mode index
x1 = real(t.x(:)); x2 = imag(t.x(:));              % target (column)
y1 = real(s.x(:))'; y2 = imag(s.x(:))';            % source (row)
ny1 = real(s.nx(:))'; ny2 = imag(s.nx(:))';        % source normal (row)
r2  = bsxfun(@minus,x1,y1).^2 + bsxfun(@minus,x2,y2).^2;
chi = 1 + r2./(2*bsxfun(@times,x1,y1));
tmp = 2./(chi+1); [K,E] = ellipke(tmp); sm = sqrt(tmp);
Qa = sm.*K;                                      % Q_{-1/2}
Qb = chi.*Qa - (chi+1).*sm.*E;                   % Q_{1/2}
for k = 1:n
  Qc = (2*k*chi.*Qb - (k-0.5).*Qa)./(k+0.5); Qa = Qb; Qb = Qc;
end
VK = Qa;                                           % Q_{n-1/2}
VE = (2*n+1).*(chi.*Qa - Qb)./(chi+1);             % = sm.*E at n=0
sqrty1x1 = sqrt(bsxfun(@rdivide,y1,x1)); x1y1 = bsxfun(@times,x1,y1);
SK = -sqrty1x1.*bsxfun(@times,ny1,x1)./x1y1;
mp = bsxfun(@minus,bsxfun(@times,ny1,x1),bsxfun(@times,ny1,y1)) ...
   + bsxfun(@minus,bsxfun(@times,ny2,x2),bsxfun(@times,ny2,y2));
SE = 2*sqrty1x1.*(bsxfun(@times,ny1,x1).*(chi-1)+mp)./r2;
A  = 1/(4*pi)*(VK.*SK + VE.*SE)*diag(s.ws);
end

function phi = inoutfun_torus(so,side)
% define indicator function for torus
my_eps = 1e-9;
if side == 'e'
  phi = @(x,y,z) ((x-so.a*cos(atan2(y,x))).^2 + (y-so.a*sin(atan2(y,x))).^2 + z.^2) > so.b^2+my_eps; 
else
  phi = @(x,y,z) ((x-so.a*cos(atan2(y,x))).^2 + (y-so.a*sin(atan2(y,x))).^2 + z.^2) < so.b^2-my_eps;
end
end

function phi = inoutfun_sphere(side)
my_eps = 1e-9;
if side == 'e'
  phi = @(x,y,z) (x.^2 + y.^2 + z.^2) > 1 + my_eps;
else
  phi = @(x,y,z) (x.^2 + y.^2 + z.^2) < 1 - my_eps;
end
end

function [f, gradf, y_src] = torus_surf_den(so, side)
switch side
  case 'e', y_src = [so.a; 0; 0]; 
  case 'i', y_src = [1.0; 1.0; 1.2]; 
end
r     = @(x) sqrt(sum((x - y_src).^2, 1));
f     = @(x) 1./(4*pi*r(x));
gradf = @(x) -(x - y_src)./(4*pi*r(x).^3);
end

function [f, gradf, y_src] = sphere_surf_den(side)
switch side
  case 'e', y_src = [0.2; 0.1; -0.3];    % inside the unit sphere
  case 'i', y_src = [1.0; 1.0; 1.2];     % outside
end
r     = @(x) sqrt(sum((x - y_src).^2, 1));
f     = @(x) 1./(4*pi*r(x));
gradf = @(x) -(x - y_src)./(4*pi*r(x).^3);
end