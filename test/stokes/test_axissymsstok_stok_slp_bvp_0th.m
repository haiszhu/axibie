clearvars; close all; format short e;
addpath('../../utils');
addpath('../../matlab');

p=16; iside=1; iclosed=0; mu=1;
lam=0.75; Z=@(t) -(1.5+cos(t)).*(-sin(lam*pi*sin(t))+1i*cos(lam*pi*sin(t)));  % c-shape (open, poles t=0,pi)
Np=6:2:36;

% 1. exact axisymmetric Stokeslet field: ring forces INSIDE the body
ysrc.x=[0.70-1.20i; 0.50-1.80i; 1.30+0.20i]; ysrc.ws=ones(size(ysrc.x));
F=[1;-0.7;0.5; 0.4;0.9;-0.6];                          % [F_rho(3); F_zeta(3)]
uexact=@(tt) StoAxiMat_even(tt,ysrc)*F;                % 2N -> [u_rho; u_zeta]

% 2. exterior meridian target grid (rho>=0), exact field on it
nx=70; gx=(1:nx)/nx*4; ny=80; gy=((1:ny)/ny*2-1)*3;
[xx,yy]=meshgrid(gx,gy); zz=xx+1i*yy;
s0=mkcurve(Z,p,Np(end));
[IN,ON]=inpolygon(real(zz),imag(zz),real(s0.x),imag(s0.x)); ii=(~IN)&(~ON);
t=[]; t.x=zz(ii(:));
ue=uexact(t); uf=nan*(1+1i)*zz; uf(ii)=ue(1:end/2)+1i*ue(end/2+1:end);

% 3. h-refinement: solve A mu = u*|_Gamma, eval A mu, error vs exact
err=nan(size(Np));
for k=1:numel(Np)
  s=mkcurve(Z,p,Np(k)); N=numel(s.x); np=s.np;
  sx=s.x(:); snx=s.nx(:); sws=s.ws(:); swxp=s.wxp(:);
  tpan=s.tpan(:); sxlo=s.Z(s.tpan(1:end-1)); sxlo=sxlo(:); sxhi=s.Z(s.tpan(2:end)); sxhi=sxhi(:);
  f=uexact(s);
  A =axss_slp_blockmat_mex(N,        sx,    p,np,sx,snx,sws,swxp,tpan,sxlo,sxhi,iside,iclosed,mu,[]);
  Ae=axss_slp_blockmat_mex(numel(t.x),t.x(:),p,np,sx,snx,sws,swxp,tpan,sxlo,sxhi,iside,iclosed,mu,[]);
  dens=A\f; tt=Ae*dens;
  u=nan*(1+1i)*zz; u(ii)=tt(1:end/2)+1i*tt(end/2+1:end);
  err(k)=max(abs(u(:)-uf(:)));
  fprintf('Stokes SLP 0th-mode Dirichlet BVP (delta mex): np=%2d, N=%d: max err vs exact = %.3e\n',Np(k),N,err(k));
end

% 4. meridian error map + convergence
figure(1); clf; semilogy(Np,err,'o-k'); xlabel('n_p'); ylabel('max err'); grid on;
title('Stokes SLP 0th-mode Dirichlet BVP: h-refinement');
figure(2); clf; imagesc(gx,gy,log10(abs(u-uf))); colorbar; hold on;
fill(real([s.x;s.x(1)]),imag([s.x;s.x(1)]),'w'); plot(s.x,'-k'); plot(ysrc.x,'r.','MarkerSize',12);
clim([-15 -8]); colormap('jet'); axis equal tight;
xlabel('\rho'); ylabel('z'); title('Stokes SLP 0th-mode: log_{10} err vs exact axi-Stokeslet');

% exportgraphics(figure(1),'axissymsstok_stok_slp0th_convergence.png','Resolution',200)
% exportgraphics(figure(2),'axissymsstok_stok_slp0th_error.png','Resolution',200)

function s=mkcurve(Z,p,np)
s.p=p; s.Z=Z; s.tpan=linspace(0,pi,np+1)'; s=quadr(s,[],'p','G');
if ~isfield(s,'np'), s.np=numel(s.tpan)-1; end
end

function G = StoAxiMat_even(t,s)
% Naive 0th-mode axisymmetric Stokes SLP, 2x2 meridian tensor (mu=1).
%   G = [G_rr G_rz; G_zr G_zz]  (2 nt x 2 ns), in-plane density stacked (sigma_rho; sigma_zeta).
%   K^{(0)}_{S,ab} = sqrt(rho'/rho)/(8 pi mu) ( SK_ab gK + SE_ab gE ) .* ws
% SK,SE the 2x2 of stokes_0th_mode_even_terms.tex (= axibie AxiStokesSpecialMat S11..S22),
% verified vs the direct 3D azimuthal integral to 1e-14.
mu = 1;
rt = real(t.x(:));    zt = imag(t.x(:));
rs = real(s.x(:)).';  zs = imag(s.x(:)).';
nt = numel(rt);       ns = numel(rs);
rho  = repmat(rt,1,ns);   rhop = repmat(rs,nt,1);
drho = rho - rhop;
dz   = repmat(zt,1,ns) - repmat(zs,nt,1);
rr2  = drho.^2 + dz.^2;
chi  = 1 + rr2./(2*rho.*rhop);   m = 2./(chi+1);
[Kc,Ec] = ellipke(m);
gK = sqrt(2./(chi+1)).*Kc;   gE = sqrt(2./(chi+1)).*Ec;
pref = (sqrt(rhop./rho)/(8*pi*mu)) .* s.ws(:).';
% 2x2 SK, SE
SKrr = (rho.^2+rhop.^2+2*dz.^2)./(rho.*rhop);   SErr = -2-4*chi+2*chi.*drho.^2./rr2;
SKrz = dz./rho;                                 SErz = 2*(drho.*dz./rr2 - dz./(2*rho));
SKzr = -dz./rhop;                               SEzr = 2*(drho.*dz./rr2 + dz./(2*rhop));
SKzz = 2*ones(size(rho));                       SEzz = 2*dz.^2./rr2;
Grr = pref.*(SKrr.*gK + SErr.*gE);
Grz = pref.*(SKrz.*gK + SErz.*gE);
Gzr = pref.*(SKzr.*gK + SEzr.*gE);
Gzz = pref.*(SKzz.*gK + SEzz.*gE);
G = [Grr Grz; Gzr Gzz];
end