clearvars; close all; format short e;
addpath('../../utils');
addpath('../../matlab');

p=16; iside=1; iclosed=0;
lam=0.75; Z=@(t) -(1.5+cos(t)).*(-sin(lam*pi*sin(t))+1i*cos(lam*pi*sin(t)));  % c-shape (open, poles t=0,pi)
Np=6:2:36;

% 1. exact axisymmetric (0th-mode) Laplace field: ring charges INSIDE the body
ysrc.x=[0.70-1.20i; 0.50-1.80i; 1.30+0.20i]; ysrc.q=[1.0; -0.7; 0.5];
uexact=@(t) lapring(t,ysrc);                           % N -> scalar potential

% 2. exterior meridian target grid (rho>=0), exact field on it
nx=70; gx=(1:nx)/nx*4; ny=80; gy=((1:ny)/ny*2-1)*3;
[xx,yy]=meshgrid(gx,gy); zz=xx+1i*yy;
s0=mkcurve(Z,p,Np(end));
[IN,ON]=inpolygon(real(zz),imag(zz),real(s0.x),imag(s0.x)); ii=(~IN)&(~ON);
t=[]; t.x=zz(ii(:));
ue=uexact(t); uf=nan*zz; uf(ii)=ue;

% 3. h-refinement: solve A sigma = u|_Gamma, eval A sigma, error vs exact
err=nan(size(Np));
for k=1:numel(Np)
  s=mkcurve(Z,p,Np(k)); N=numel(s.x); np=s.np;
  sx=s.x(:); snx=s.nx(:); sws=s.ws(:); swxp=s.wxp(:);
  tpan=s.tpan(:); sxlo=s.Z(s.tpan(1:end-1)); sxlo=sxlo(:); sxhi=s.Z(s.tpan(2:end)); sxhi=sxhi(:);
  f=uexact(s);
  A =axls_slp_blockmat_mex(N,         sx,    p,np,sx,snx,sws,swxp,tpan,sxlo,sxhi,iside,iclosed,[]);
  Ae=axls_slp_blockmat_mex(numel(t.x),t.x(:),p,np,sx,snx,sws,swxp,tpan,sxlo,sxhi,iside,iclosed,[]);
  dens=A\f; tt=Ae*dens;
  u=nan*zz; u(ii)=tt;
  err(k)=max(abs(u(:)-uf(:)));
  fprintf('Laplace SLP 0th-mode Dirichlet BVP (delta mex): np=%2d, N=%d: max err vs exact = %.3e\n',Np(k),N,err(k));
end

% 4. meridian error map + convergence
figure(1); clf; semilogy(Np,err,'o-k'); xlabel('n_p'); ylabel('max err'); grid on;
title('Laplace SLP 0th-mode Dirichlet BVP: h-refinement');
figure(2); clf; imagesc(gx,gy,log10(abs(u-uf))); colorbar; hold on;
fill(real([s.x;s.x(1)]),imag([s.x;s.x(1)]),'w'); plot(s.x,'-k'); plot(ysrc.x,'r.','MarkerSize',12);
clim([-15 -11]); colormap('jet'); axis equal tight;
xlabel('\rho'); ylabel('z'); title('Laplace SLP 0th-mode: log_{10} err vs exact axi-ring potential');

% cmp = getPyPlot_cMap('rainbow', [], [], '"/Users/hzhu/.pyenv/versions/3.11.13/bin/python"');
% colormap(cmp)

% exportgraphics(figure(1),'axissymslap_lap_slp0th_convergence.png','Resolution',200)
% exportgraphics(figure(2),'axissymslap_lap_slp0th_error.png','Resolution',200)

function s=mkcurve(Z,p,np)
s.p=p; s.Z=Z; s.tpan=linspace(0,pi,np+1)'; s=quadr(s,[],'p','G');
if ~isfield(s,'np'), s.np=numel(s.tpan)-1; end
end

function u = lapring(t,ysrc)
% axisymmetric (0th-mode) Laplace potential of ring charges ysrc (meridian rho+iz, charge q), via the
% elliptic-K ring formula:  Phi(rho,z) = sum_j q_j/(2 pi^2) K(k_j^2)/sqrt((rho+a_j)^2+(z-b_j)^2),
% with k_j^2 = 4 rho a_j / ((rho+a_j)^2+(z-b_j)^2).
rho=real(t.x(:)); z=imag(t.x(:)); u=zeros(numel(rho),1);
for j=1:numel(ysrc.x)
  a=real(ysrc.x(j)); b=imag(ysrc.x(j));
  d2=(rho+a).^2+(z-b).^2; k2=4*rho*a./d2;
  u=u+ysrc.q(j)/(2*pi^2)*ellipke(k2)./sqrt(d2);
end
end
