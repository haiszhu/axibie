addpath('../../matlab/')
addpath('../../utils/')
addpath('../../external/fmm3d/matlab/')

clear all

% set parameters
p    = 16;
np   = 8; % 1.42e-06
np   = 12; % 1.1e-08
np   = 16; % 3.78e-10
% np   = 20; % 1.16e-09, it is unclear why... need to figure out
M    = 4*np;
nphi = 2*M+1;
type = 'torus';
side = 'e';
N    = np*p;
k    = 5;
iside = double(side=='e');

% generate torus/sphere geometry (full 3d + axisymmetric)
so.a = 1.0; so.b = 0.5;
s = [];
s.p = p; s.Z  = @(t) (so.a+so.b*cos(t-pi/2)) + 1i*(so.b*sin(t-pi/2)); s.tpan = linspace(0, 2*pi, np+1)'; 
s   = quadr(s,[],'p','G');
s3d = axisym_to_3d_quadrature( real(s.x), imag(s.x), s.ws, real(s.nx), imag(s.nx), nphi);

% define density function for GRF from a source
[u_vel, f_trac]= torus_surf_den(so,side);

% all on surface node
sx = s3d.x; snx = s3d.nx; sw = s3d.w;

sigma = f_trac(sx,snx); tau = -u_vel(sx); % densities for GRF 

% exterior GRF test
nt = 50000;
rng = 4;
t3d.x = [0.8520;0.4919;0.4819] + 2.0*(rand(3,nt)-1/2); % random targets, shrink 2*sqrt(sum(pan{k}.w)), or increase nt to observe closer targets...
phi = inoutfun(so,side);
inout_flag = phi(t3d.x(1,:),t3d.x(2,:),t3d.x(3,:));
t3d_inout.x = t3d.x(:,inout_flag);
X = reshape(s3d.x(1,:), numel(s.x), nphi); Y = reshape(s3d.x(2,:), numel(s.x), nphi); Z = reshape(s3d.x(3,:), numel(s.x), nphi); X=X([1:end 1],[1:end 1]); Y=Y([1:end 1],[1:end 1]); Z=Z([1:end 1],[1:end 1]);

% solution
fhom = u_vel(t3d_inout.x);

% naive evaluation
eps = 1e-14;
slpval = Sto3dSLPfmm_il(t3d_inout,s3d,sigma,eps);   % (3*Mi,1) interleaved column
dlpval = Sto3dDLPfmm_il(t3d_inout,s3d,tau,eps);
u3d_inout = reshape(slpval + dlpval, 3, []);        % -> 3 x Mi, matches fhom

% reorder FFT output to match
modes = -M:M;
th = atan2(sx(2,:), sx(1,:));
c2y = @(v)[ v(1,:).*cos(th)+v(2,:).*sin(th);           % vρ
           -v(1,:).*sin(th)+v(2,:).*cos(th);           % vθ   <- the 3 scalars
            v(3,:) ];                                  % vz
sig_cyl = c2y(sigma);  tau_cyl = c2y(tau);
sigma_hat = zeros(2*M+1,N,3); tau_hat = zeros(2*M+1,N,3);
for c = 1:3                                            % Laplace block, verbatim, per component
  sigma_hat(:,:,c) = fftshift(fft(reshape(sig_cyl(c,:),N,nphi)',2*M+1,1)/nphi,1);
  tau_hat(:,:,c)   = fftshift(fft(reshape(tau_cyl(c,:),N,nphi)',2*M+1,1)/nphi,1);
end

% error visualization (per-target velocity-vector norm)
if side == 'i', err3d = vecnorm(u3d_inout - fhom, 2, 1); end   % source outside -> GRF reproduces +u_ex
if side == 'e', err3d = vecnorm(u3d_inout + fhom, 2, 1); end   % source inside  -> GRF reproduces -u_ex
fprintf('N=%d (num of source), naive GRF test at %d target pt: u err = %.3g\n', size(s3d.x,2), numel(t3d_inout.x)/3, max(err3d));
figure(1),clf,
subplot(1,3,1)
scatter3(t3d_inout.x(1,:),t3d_inout.x(2,:),t3d_inout.x(3,:),7,log10(err3d),'filled'),colorbar,
hold on, h=surf(X,Y,Z); set(h,'ambientstrength',0.7,'facelighting','gouraud'); shading interp; set(h,'FaceColor','w'); lightangle(45,0); axis equal vis3d; set(gca,'clipping','off'); xlabel('x'); ylabel('y'); zlabel('z'); colorbar;
caxis([-16 2]),title(sprintf('naive full-3D (max err = %.1e)', max(err3d)))
view(150,55)
if ismac
  pyCmd = '"/Users/hzhu/.pyenv/versions/3.11.13/bin/python"';
  % pyCmd = '/Users/hzhu/.venvs/mpl/bin/python';   % macOS venv alt
else
  pyCmd = ['"' fullfile(fileparts(mfilename('fullpath')),'..','..','..','.venv','bin','python') '"'];   % repo uv venv (matplotlib)
end
cmp = getPyPlot_cMap('rainbow', [], [], pyCmd); colormap(cmp)
colorbar off; ax=gca; cb=colorbar; drawnow; ax.Position=ax.Position; p=cb.Position; cb.Position=[p(1) p(2)+p(4)/4 p(3) p(4)/2];

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
thk = tk_inout.theta(:).';
mu = 1;
uk_skip_k = zeros(3, nk);
As = axss_slp_blockmat_nmode_mex(nk, tk_inout.x(:).', s.p, s.np, s.x, s.nx, s.ws, s.wxp, s.tpan, s.xlo, s.xhi, M, iside, 1, mu, [], 1, k);
Ad = axss_dlp_blockmat_nmode_mex(nk, tk_inout.x(:).', s.p, s.np, s.x, s.nx, s.ws, s.wxp, s.tpan, s.xlo, s.xhi, M, iside, 1, mu, [], 1, k);
for jj = 1:numel(modes)
  n = modes(jj);
  Sblk = As(:,:,abs(n)+1);  Dblk = Ad(:,:,abs(n)+1);      % mex stores modes 0..M
  if n<0, Sblk = conj(Sblk); Dblk = conj(Dblk); end       % A(-n)=conj(A(n)) (Stokes kernel not even in n)
  vk = reshape(Sblk*sigma_hat(jj,:).' + Dblk*tau_hat(jj,:).', nk, 3).';   % 3 x nk cyl (block-ordered out)
  uk_skip_k = uk_skip_k + vk.*exp(1i*n*tk_inout.theta(:).');
end
uk_skip_k = real(uk_skip_k);

% add back naive k-th panel contribution
uk_only_k = zeros(3, nk);
idxk = (k-1)*s.p + (1:s.p);
sk = []; sk.x = s.x(idxk); sk.nx = s.nx(idxk); sk.ws = s.ws(idxk);
for jj = 1:numel(modes)
  n  = modes(jj);
  Snk = StoSLPAxiMattmp(tk_inout, sk, n);
  Dnk = StoDLPAxiMattmp(tk_inout, sk, n);
  sigk = reshape(squeeze(sigma_hat(jj,idxk,:)).', [], 1);
  tauk = reshape(squeeze(tau_hat(jj,idxk,:)).',  [], 1);
  vntmp_k = reshape(Snk*sigk + Dnk*tauk, 3, nk);
  uk_only_k = uk_only_k + vntmp_k.*exp(1i*n*tk_inout.theta(:).');
end
uk_only_k = real(uk_only_k);

% combine + cyl -> cart, then compare (mirrors your stage-5 tail)
uk2 = uk_skip_k + uk_only_k;                              % 3 x nk cylindrical
uk2 = [uk2(1,:).*cos(thk)-uk2(2,:).*sin(thk); uk2(1,:).*sin(thk)+uk2(2,:).*cos(thk); uk2(3,:)];
fhomk = fhom(:,sqn_idx);
if side=='i', errk2 = vecnorm(uk2-fhomk,2,1); else, errk2 = vecnorm(uk2+fhomk,2,1); end
fprintf('without k-th panel close-eval corrected: max err at %d close targets = %.3g\n', nk, max(errk2));
figure(1),
subplot(1,3,2)
scatter3(tk3d_inout.x(1,:), tk3d_inout.x(2,:), tk3d_inout.x(3,:), 12, log10(errk2), 'filled'); colorbar, colormap(cmp)
hold on, h=surf(X,Y,Z); set(h,'ambientstrength',0.7,'facelighting','gouraud'); shading interp; set(h,'FaceColor','w'); lightangle(45,0); axis equal vis3d; set(gca,'clipping','off'); xlabel('x'); ylabel('y'); zlabel('z'); colorbar;
axis equal vis3d; set(gca,'clipping','off'); xlabel('x'); ylabel('y'); zlabel('z'); caxis([-16 2]);
title(sprintf('naive eval (panel k=%d / %d)', k, np))
view(150,55)
colorbar off; ax=gca; cb=colorbar; drawnow; ax.Position=ax.Position; p=cb.Position; cb.Position=[p(1) p(2)+p(4)/4 p(3) p(4)/2];

% use close k-th panel contribution instead
skx = s.x((k-1)*s.p+(1:s.p));
[As3d, Ad3d] = axps_closestoksdlp_panel_mex(nk, tk_inout.x, tk3d_inout.x, [], s.p, nphi, skx, [], [], [], ...
                                          s.xlo(k), s.xhi(k), [], [], [], [], [], M, iside, [], mu, ...
                                          [], [], [], [], [], []);
cols = (k-1)*s.p;
ring3 = reshape(3*(reshape((cols+(1:s.p)).' + (0:nphi-1)*N,1,[])-1) + (1:3).', [], 1);
sig3  = sigma(:);  tau3 = tau(:);
uk_only_k_close = reshape(As3d*sig3(ring3) + Ad3d*tau3(ring3), 3, nk);
uk_skip_k_cart = [uk_skip_k(1,:).*cos(thk) - uk_skip_k(2,:).*sin(thk); ...
                  uk_skip_k(1,:).*sin(thk) + uk_skip_k(2,:).*cos(thk); ...
                  uk_skip_k(3,:)];
uk3 = uk_skip_k_cart + uk_only_k_close;
if side=='i', errk3 = vecnorm(uk3-fhomk,2,1); else, errk3 = vecnorm(uk3+fhomk,2,1); end
fprintf('with k-th panel close-eval (panel wrapper): max err at %d close targets = %.3g\n', nk, max(errk3));
figure(1),
subplot(1,3,3)
scatter3(tk3d_inout.x(1,:), tk3d_inout.x(2,:), tk3d_inout.x(3,:), 12, log10(errk3), 'filled'); colorbar, colormap(cmp)
hold on, h=surf(X,Y,Z); set(h,'ambientstrength',0.7,'facelighting','gouraud'); shading interp; set(h,'FaceColor','w'); lightangle(45,0); axis equal vis3d; set(gca,'clipping','off'); xlabel('x'); ylabel('y'); zlabel('z'); colorbar;
axis equal vis3d; set(gca,'clipping','off'); xlabel('x'); ylabel('y'); zlabel('z'); caxis([-16 2]);
title(sprintf('with close eval (panel k=%d / %d)', k, np))
view(150,55)
colorbar off; ax=gca; cb=colorbar; drawnow; ax.Position=ax.Position; p=cb.Position; cb.Position=[p(1) p(2)+p(4)/4 p(3) p(4)/2];

% keyboard
% exportgraphics(figure(1),'axissymsstok_stok_greens_identity.png','Resolution',200)

function A = StoSLPAxiMattmp(t,s,n)
rho  = real(t.x(:));   zt = imag(t.x(:));           % target (column)
rhop = real(s.x(:)).'; zs = imag(s.x(:)).';         % source (row)
zh   = zt - zs;                                      % nt x ns
chi  = 1 + ((rho-rhop).^2 + zh.^2)./(2*rho.*rhop);
nn = abs(n);                                         % carriers use |n|; the 9 coeffs keep signed n
tmp  = 2./(chi+1); [K,E] = ellipke(tmp); sm = sqrt(tmp);
Qa = sm.*K; Qb = chi.*Qa - (chi+1).*sm.*E;         % Q_{-1/2}, Q_{1/2}
for k = 1:nn
  Qc = (2*k*chi.*Qb - (k-0.5).*Qa)./(k+0.5); Qa = Qb; Qb = Qc;
end
VK = Qa; VE = (2*nn+1).*(chi.*Qa - Qb)./(chi+1);    % Q_{|n|-1/2}, V_E
ip = 1/pi; rr = rho.*rhop; rrt = (rho.*rhop).^(-3/2); srt = 1./sqrt(rho.*rhop);
n2 = n^2; fn2 = 4*n2-1; cm1 = chi-1; cp1 = chi+1;
SK = cell(1,9); SE = cell(1,9);
SK{1} = -(1/8).*rhop.*ip.*rrt.*(rho.^2+rhop.^2-4*chi.*rr);                                % rr
SE{1} = -(1/8).*rhop.*ip.*rrt.*(1./cm1).*(1./fn2).*( 2*rr+4*n2.*rr+chi.*rho.^2+chi.*rhop.^2 ...
        -4*chi.^2.*rr+4*n2.*chi.^2.*rr-4*chi.*n2.*rho.^2-4*chi.*n2.*rhop.^2 );
SK{2} = -(1/4).*n.*srt.*ip.*1i.*(rho-chi.*rhop);                                          % rphi
SE{2} = -(3/4).*n.*rhop.*srt.*ip.*1i.*cp1.*(1./fn2);
SK{3} =  (1/8).*rhop.^2.*ip.*rrt.*zh;                                                     % rz
SE{3} =  (1/8).*rhop.*ip.*rrt.*(1./cm1).*zh.*(rho-chi.*rhop);
SK{4} =  (1/4).*n.*rhop.^2.*ip.*rrt.*1i.*(rhop-chi.*rho);                                 % phir
SE{4} =  (3/4).*n.*rhop.*srt.*ip.*1i.*cp1.*(1./fn2);
SK{5} =  (1/2).*chi.*rhop.*srt.*ip;                                                       % phiphi
SE{5} = -(1/2).*rhop.*srt.*ip.*cp1.*(n2-1).*(1./fn2);
SK{6} = -(1/4).*1i.*n.*rhop.^2.*ip.*rrt.*zh;                                              % phiz
SE{6} =  zeros(size(chi));
SK{7} = -(1/8).*srt.*ip.*zh;                                                              % zr
SE{7} = -(1/8).*rhop.*ip.*rrt.*(1./cm1).*zh.*(rhop-chi.*rho);
SK{8} = -(1/4).*1i.*n.*srt.*ip.*zh;                                                       % zphi
SE{8} =  zeros(size(chi));
SK{9} =  (1/4).*rhop.*srt.*ip;                                                            % zz
SE{9} =  (1/8).*rhop.*ip.*rrt.*(1./cm1).*zh.^2;
ws = s.ws(:).';
G = cell(1,9); for e = 1:9, G{e} = (SK{e}.*VK + SE{e}.*VE).*ws; end
[nt,ns] = size(G{1}); A = zeros(3*nt,3*ns);            % node-interleaved (xyz,xyz,...) rows & cols
A(1:3:end,1:3:end)=G{1}; A(1:3:end,2:3:end)=G{2}; A(1:3:end,3:3:end)=G{3};
A(2:3:end,1:3:end)=G{4}; A(2:3:end,2:3:end)=G{5}; A(2:3:end,3:3:end)=G{6};
A(3:3:end,1:3:end)=G{7}; A(3:3:end,2:3:end)=G{8}; A(3:3:end,3:3:end)=G{9};
end

function A = StoDLPAxiMattmp(t,s,n)
rho  = real(t.x(:));    zt = imag(t.x(:));           % target (column)
rhop = real(s.x(:)).';  zs = imag(s.x(:)).';         % source (row)
nrp  = real(s.nx(:)).'; nzp = imag(s.nx(:)).';       % source normal (row)
zh   = zt - zs;                                       % nt x ns
chi  = 1 + ((rho-rhop).^2 + zh.^2)./(2*rho.*rhop);
nn = abs(n);                                          % carriers use |n|; the 9 coeffs keep signed n
tmp  = 2./(chi+1); [K,E] = ellipke(tmp); sm = sqrt(tmp);
Qa = sm.*K; Qb = chi.*Qa - (chi+1).*sm.*E;          % Q_{-1/2}, Q_{1/2}
for k = 1:nn
  Qc = (2*k*chi.*Qb - (k-0.5).*Qa)./(k+0.5); Qa = Qb; Qb = Qc;
end
VK = Qa; VE = (2*nn+1).*(chi.*Qa - Qb)./(chi+1);      % Q_{|n|-1/2}, V_E
% --- inline dlp9coef(n,chi,rho,rhop,zh,nrp,nzp): shared temps + the (S_K,S_E) of each of 9 entries ---
t2=nrp.*rho; t3=nrp.*rhop; t4=rho.*rhop; t5=nzp.*zh; t6=chi+1.0; t7=chi.^2; t8=chi.^3;
t10=n.^2; t11=rho.^2; t13=rhop.^2; t15=zh.^2; t17=1.0./pi; t20=chi-1.0; t21=-rhop;
t9=t7.^2; t16=t2.*3.0; t18=rho.*t5; t19=rhop.*t5; t22=t10.*4.0; t23=-t3; t24=-t5;
t25=rho.*t2; t26=rhop.*t3; t27=t3.*t13; t28=chi.*t2.*4.0; t29=rhop.*t2.*4.0; t30=chi.*rhop.*t2.*2.0;
t31=t5.*t13; t32=1.0./t6; t33=t7-1.0; t34=t4.*t5.*4.0; t35=t2.*1i; t36=t3.*1i; t37=t3.*2.0i;
t38=t5.*1i; t42=rho+t21; t47=t2.*t13.*7.0; t48=1.0./t20; t53=chi.*t2.*2.0i; t58=rhop.*t2.*t7.*2.0;
t62=1.0./t4.^(5.0./2.0); t67=rhop.*t2.*t10.*4.0i; t39=rho.*t16; t40=t11.*t16; t41=t26.*3.0;
t43=-t29; t44=-t18; t45=-t19; t46=t4.*t16; t49=t48.^2; t50=-t34; t51=-t36; t52=chi.*t35;
t54=t19.*1i; t55=t42.^2; t56=rhop.*t53; t57=t13.*t24; t59=t22-1.0; t60=t25.*1i; t61=t26.*1i;
t64=1.0./t42; t65=1.0./t33;
t68=t10.*t19.*4.0i; t69=t10.*t26.*4.0i; t70=1.0./t59; t73=t38+t51+t52;
SK = cell(1,9); SE = cell(1,9);
SK{1} = rhop.*t17.*t62.*t65.*(t40+t47+t50+chi.*t31-chi.*t2.*t4.*1.1e+1+chi.*t5.*t11+chi.*t13.*t23-chi.*t10.*t31.*4.0+chi.*t22.*t27+t2.*t4.*t8.*8.0+t4.*t5.*t7.*2.0-t2.*t7.*t11.*2.0-t2.*t7.*t13.*4.0-t2.*t10.*t13.*4.0+t4.*t5.*t22+chi.*t2.*t4.*t10.*8.0-chi.*t5.*t10.*t11.*4.0-t2.*t7.*t10.*t11.*4.0-t2.*t7.*t10.*t13.*8.0+t2.*t4.*t8.*t22+t4.*t5.*t7.*t22).*(-1.0./8.0);                                                                                            % rr S_K
SE{1} = (rhop.*t17.*t32.*t49.*t62.*t70.*(t27.*3.0-t31.*3.0+t46-t5.*t11.*3.0+t7.*t27-t10.*t27.*1.2e+1+t10.*t31.*1.2e+1+t7.*t57+chi.*t4.*t5.*1.0e+1-chi.*t2.*t11.*6.0-chi.*t2.*t13.*1.6e+1+t2.*t4.*t7.*1.7e+1-t2.*t4.*t9.*8.0-t2.*t4.*t10.*2.4e+1-t4.*t5.*t8.*2.0+t2.*t8.*t11.*2.0+t2.*t8.*t13.*4.0+t5.*t10.*t11.*1.2e+1+t7.*t11.*t24-t7.*t10.*t27.*4.0+t7.*t22.*t31-chi.*t4.*t5.*t10.*4.0e+1+chi.*t2.*t10.*t11.*2.4e+1+chi.*t2.*t10.*t13.*6.4e+1-t2.*t4.*t7.*t10.*4.4e+1+t2.*t4.*t9.*t10.*2.0e+1+t4.*t5.*t8.*t10.*8.0-t2.*t8.*t10.*t11.*8.0-t2.*t8.*t10.*t13.*1.6e+1+t5.*t7.*t11.*t22))./8.0;   % rr S_E
SK{2} = n.*1.0./t4.^(3.0./2.0).*t17.*(t54+t56-t60-t61).*(-3.0./4.0);                                    % rphi S_K
SE{2} = (n.*1.0./t4.^(3.0./2.0).*t17.*t48.*t70.*(t18.*1i+t67-chi.*t19.*1i+chi.*t60+chi.*t61+chi.*t68-rhop.*t2.*4.0i-t10.*t18.*4.0i+t7.*t67-chi.*t10.*t25.*4.0i-chi.*t10.*t26.*4.0i+rhop.*t2.*t7.*2.0i))./4.0;   % rphi S_E
SK{3} = rhop.*t17.*t62.*t65.*zh.*(t18+t43+t58+chi.*t25+chi.*t26+chi.*t45-t10.*t18.*4.0-chi.*t10.*t25.*4.0-chi.*t10.*t26.*4.0+chi.*t19.*t22+rhop.*t2.*t22+rhop.*t2.*t7.*t22).*(-1.0./8.0);                        % rz S_K
SE{3} = (rhop.*t17.*t32.*t49.*t62.*zh.*(t19.*-3.0+t39+t41+chi.*t18.*4.0+t7.*t25+t7.*t26+t7.*t45-chi.*rhop.*t2.*1.0e+1+rhop.*t2.*t8.*2.0))./8.0;   % rz S_E
SK{4} = n.*rhop.*1.0./t4.^(3.0./2.0).*t17.*(-t37+t38+t53).*(3.0./4.0);                                  % phir S_K
SE{4} = n.*t13.*t17.*t48.*t62.*t70.*(t25.*-3.0i+t54+t56-t61-t68+t69-chi.*t18.*1i+t7.*t25.*2.0i+chi.*t10.*t18.*4.0i+t7.*t10.*t25.*4.0i-chi.*rhop.*t2.*t10.*8.0i).*(-1.0./4.0);   % phir S_E
SK{5} = rhop.*1.0./t4.^(3.0./2.0).*t17.*(t5+t23+t28-t3.*t10.*2.0+t5.*t10.*2.0+chi.*t2.*t10.*2.0).*(-1.0./4.0);   % phiphi S_K
SE{5} = (rhop.*1.0./t4.^(3.0./2.0).*t17.*t48.*t70.*(t16+chi.*t3+chi.*t24-t2.*t7.*4.0-t2.*t10.*6.0-chi.*t3.*t10.*4.0+chi.*t5.*t22+t2.*t7.*t10.*1.0e+1))./4.0;   % phiphi S_E
SK{6} = n.*t2.*t13.*t17.*t62.*zh.*7.5e-1i;                                                              % phiz S_K
SE{6} = n.*t13.*t17.*t48.*t62.*t73.*zh.*(-1.0./4.0);                                                    % phiz S_E
SK{7} = (rhop.*t17.*t62.*t65.*zh.*(t19-t25.*3.0+t30+chi.*t44+t3.*t21-t10.*t19.*4.0+t7.*t25.*2.0+t22.*t26+chi.*t18.*t22+t7.*t22.*t25-chi.*rhop.*t2.*t10.*8.0))./8.0;   % zr S_K
SE{7} = (rhop.*t17.*t32.*t49.*t62.*zh.*(t18.*3.0-t58-chi.*t19.*4.0+chi.*t25.*6.0+chi.*t26.*4.0-rhop.*t2.*6.0+t7.*t18-t8.*t25.*2.0))./8.0;   % zr S_E
SK{8} = n.*t2.*1.0./t4.^(3.0./2.0).*t17.*zh.*7.5e-1i;                                                   % zphi S_K
SE{8} = n.*1.0./t4.^(3.0./2.0).*t17.*t48.*t73.*zh.*(-1.0./4.0);                                         % zphi S_E
SK{9} = (rhop.*t15.*t17.*t59.*t62.*t65.*(t5+t23+chi.*t2))./8.0;                                         % zz S_K
SE{9} = (rhop.*t15.*t17.*t32.*t49.*t62.*(t16-chi.*t3.*4.0+chi.*t5.*4.0+t2.*t7))./8.0;                   % zz S_E
ws = s.ws(:).';
G = cell(1,9); for e = 1:9, G{e} = (SK{e}.*VK + SE{e}.*VE).*ws; end
[nt,ns] = size(G{1}); A = zeros(3*nt,3*ns);            % node-interleaved (xyz,xyz,...) rows & cols
A(1:3:end,1:3:end)=G{1}; A(1:3:end,2:3:end)=G{2}; A(1:3:end,3:3:end)=G{3};
A(2:3:end,1:3:end)=G{4}; A(2:3:end,2:3:end)=G{5}; A(2:3:end,3:3:end)=G{6};
A(3:3:end,1:3:end)=G{7}; A(3:3:end,2:3:end)=G{8}; A(3:3:end,3:3:end)=G{9};
end

function phi = inoutfun(so,side)
% define indicator function for torus
my_eps = 1e-9;
if side == 'e'
  phi = @(x,y,z) ((x-so.a*cos(atan2(y,x))).^2 + (y-so.a*sin(atan2(y,x))).^2 + z.^2) > so.b^2+my_eps; 
else
  phi = @(x,y,z) ((x-so.a*cos(atan2(y,x))).^2 + (y-so.a*sin(atan2(y,x))).^2 + z.^2) < so.b^2-my_eps;
end
end

function [u, f, y_force] = torus_surf_den(so, side)
rng(4); pt_force = rand(3,1); pt_force = 10*pt_force/norm(pt_force); 
switch side
  case 'e'
    y_force.x = [so.a; 0; 0]; 
  case 'i' 
    y_force.x = [1.0; 1.0; 1.2];
end
y_force.w = 1;
u = @(x)    reshape(Sto3dSLPmat_il(struct('x',x),y_force)*pt_force,3,[]);
f = @(x,nx) reshape(Sto3dSLPnmat_il(struct('x',x,'nx',nx),y_force)*pt_force,3,[]);

function T = Sto3dSLPnmat_il(t,s)
[~, T] = Sto3dSLPmat_il(t,s);
end

end