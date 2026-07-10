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

% 3. h-refinement: solve (D_+ + S) tau = u*|_Gamma (D_+ = exterior limit, jump included), eval, error
err=nan(size(Np));
for k=1:numel(Np)
  s=mkcurve(Z,p,Np(k)); N=numel(s.x); np=s.np;
  sx=s.x(:); snx=s.nx(:); sws=s.ws(:); swxp=s.wxp(:);
  tpan=s.tpan(:); sxlo=s.Z(s.tpan(1:end-1)); sxlo=sxlo(:); sxhi=s.Z(s.tpan(2:end)); sxhi=sxhi(:);
  f=uexact(s);
  im=[1:N,2*N+1:3*N];                                  % meridian (rho,z) rows/cols of the mode-0 3x3 stack
  % FROZEN single-mode path (uncomment to compare): A =axss_dlp_blockmat_mex(N,        sx,    p,np,sx,snx,sws,swxp,tpan,sxlo,sxhi,iside,iclosed,mu,[]) ...
  %   +axss_slp_blockmat_mex(N,        sx,    p,np,sx,snx,sws,swxp,tpan,sxlo,sxhi,iside,iclosed,mu,[]);
  A3=axp_modemat_setup_mex(2,3,mu,1,1,p,np,0,iside,iclosed,[1;N+1],[1;np+2],N,np+1, ...
      sx,snx,sws,swxp,tpan,0,[1;1],zeros(3,0),zeros(3,0),eye(3),zeros(3,1),3*N,3*N) ...
    + axp_modemat_setup_mex(2,1,mu,1,1,p,np,0,iside,iclosed,[1;N+1],[1;np+2],N,np+1, ...
      sx,snx,sws,swxp,tpan,0,[1;1],zeros(3,0),zeros(3,0),eye(3),zeros(3,1),3*N,3*N);
  A=A3(im,im);
  Mt=numel(t.x); T3=[real(t.x(:)).'; zeros(1,Mt); imag(t.x(:)).'];   % meridian targets in 3D (phi=0)
  % FROZEN single-mode path (uncomment to compare): Ae=axss_dlp_blockmat_mex(numel(t.x),t.x(:),p,np,sx,snx,sws,swxp,tpan,sxlo,sxhi,iside,iclosed,mu,[]) ...
  %   +axss_slp_blockmat_mex(numel(t.x),t.x(:),p,np,sx,snx,sws,swxp,tpan,sxlo,sxhi,iside,iclosed,mu,[]);
  Ae3=axp_modemat_setup_mex(2,3,mu,3,1,p,np,0,iside,iclosed,[1;N+1],[1;np+2],N,np+1, ...
      sx,snx,sws,swxp,tpan,Mt,[1;Mt+1],T3,zeros(3,Mt),eye(3),zeros(3,1),3*Mt,3*N) ...
    + axp_modemat_setup_mex(2,1,mu,3,1,p,np,0,iside,iclosed,[1;N+1],[1;np+2],N,np+1, ...
      sx,snx,sws,swxp,tpan,Mt,[1;Mt+1],T3,zeros(3,Mt),eye(3),zeros(3,1),3*Mt,3*N);
  Ae=Ae3([1:Mt,2*Mt+1:3*Mt],im);
  tau=A\f; tt=Ae*tau;
  u=nan*(1+1i)*zz; u(ii)=tt(1:end/2)+1i*tt(end/2+1:end);
  err(k)=max(abs(u(:)-uf(:)));
  fprintf('Stokes combined (D_+ + S) 0th-mode Dirichlet BVP (delta mex): np=%2d, N=%d: max err vs exact = %.3e\n',Np(k),N,err(k));
end

% 4. meridian error map + convergence
figure(1); clf; semilogy(Np,err,'o-k'); xlabel('n_p'); ylabel('max err'); grid on;
title('Stokes combined (D_+ + S) 0th-mode Dirichlet BVP (delta mex): h-refinement');
figure(2); clf; imagesc(gx,gy,log10(abs(u-uf))); colorbar; hold on;
fill(real([s.x;s.x(1)]),imag([s.x;s.x(1)]),'w'); plot(s.x,'-k'); plot(ysrc.x,'r.','MarkerSize',12);
clim([-15 -8]); colormap('jet'); axis equal tight;
xlabel('\rho'); ylabel('z'); title('Stokes combined (D_+ + S) 0th-mode: log_{10} err vs exact axi-Stokeslet');

% exportgraphics(figure(1),'axissymsstok_stok_dlp0th_convergence.png','Resolution',200)
% exportgraphics(figure(2),'axissymsstok_stok_dlp0th_error.png','Resolution',200)

function s=mkcurve(Z,p,np)
s.p=p; s.Z=Z; s.tpan=linspace(0,pi,np+1)'; s=quadr(s,[],'p','G');
if ~isfield(s,'np'), s.np=numel(s.tpan)-1; end
end
