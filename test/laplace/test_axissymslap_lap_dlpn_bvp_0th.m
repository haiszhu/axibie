clearvars; close all; format short e;
addpath('/Users/hzhu/Documents/Github/axibie/utils');
addpath('/Users/hzhu/Documents/Github/AxiStokes3D/utils');
addpath('/Users/hzhu/Documents/Github/AxiStokes3D/matlab');

p=16; iside=1; iclosed=0;
lam=0.75; Z=@(t) -(1.5+cos(t)).*(-sin(lam*pi*sin(t))+1i*cos(lam*pi*sin(t)));  % c-shape (open, poles t=0,pi)
Np=6:2:36;

% 1. exact axisymmetric (0th-mode) Laplace field: ring charges INSIDE the body
ysrc.x=[0.70-1.20i; 0.50-1.80i; 1.30+0.20i]; ysrc.q=[1.0; -0.7; 0.5];
uexact=@(t) lapring(t,ysrc);                            % N -> scalar potential
gexact=@(t) lapringn(t,ysrc);                           % N -> dn u (target needs t.nx)

% 2. exterior meridian target grid (rho>=0), exact field on it
nx=70; gx=(1:nx)/nx*4; ny=80; gy=((1:ny)/ny*2-1)*3;
[xx,yy]=meshgrid(gx,gy); zz=xx+1i*yy;
s0=mkcurve(Z,p,Np(end));
[IN,ON]=inpolygon(real(zz),imag(zz),real(s0.x),imag(s0.x)); ii=(~IN)&(~ON);
t=[]; t.x=zz(ii(:));
ue=uexact(t); uf=nan*zz; uf(ii)=ue;

% 3. h-refinement: combined-field Neumann (-1/2 I + S' + D') sigma = dn u, eval (S+D)[sigma]
err=nan(size(Np));
for k=1:numel(Np)
  s=mkcurve(Z,p,Np(k)); N=numel(s.x); np=s.np;
  sx=s.x(:); snx=s.nx(:); sws=s.ws(:); swxp=s.wxp(:);
  tpan=s.tpan(:); sxlo=s.Z(s.tpan(1:end-1)); sxlo=sxlo(:); sxhi=s.Z(s.tpan(2:end)); sxhi=sxhi(:);
  g=gexact(s);
  Sp=axls_slpn_blockmat_mex(N,sx,snx,p,np,sx,snx,sws,swxp,tpan,sxlo,sxhi,iside,iclosed,[]);
  Dp=axls_dlpn_blockmat_mex(N,sx,snx,p,np,sx,snx,sws,swxp,tpan,sxlo,sxhi,iside,iclosed,[]);
  dens=(Sp+Dp)\g;
  Se=axls_slp_blockmat_mex(numel(t.x),t.x(:),p,np,sx,snx,sws,swxp,tpan,sxlo,sxhi,iside,iclosed,[]);
  De=axls_dlp_blockmat_mex(numel(t.x),t.x(:),p,np,sx,snx,sws,swxp,tpan,sxlo,sxhi,iside,iclosed,[]);
  tt=(Se+De)*dens;
  u=nan*zz; u(ii)=tt;
  err(k)=max(abs(u(:)-uf(:)));
  fprintf('Laplace combined (S''+D'') 0th-mode Neumann BVP (delta mex): np=%2d, N=%d: max err vs exact = %.3e\n',Np(k),N,err(k));
end

% 4. meridian error map + convergence
figure(1); clf; semilogy(Np,err,'o-k'); xlabel('n_p'); ylabel('max err'); grid on;
title('Laplace combined (S''+D'') 0th-mode Neumann BVP: h-refinement');
figure(2); clf; imagesc(gx,gy,log10(abs(u-uf))); colorbar; hold on;
fill(real([s.x;s.x(1)]),imag([s.x;s.x(1)]),'w'); plot(s.x,'-k'); plot(ysrc.x,'r.','MarkerSize',12);
clim([-16 -11]); colormap('jet'); axis equal tight;
xlabel('\rho'); ylabel('z'); title('Laplace combined (S''+D'') 0th-mode: log_{10} err vs exact axi-ring potential');

% cmp = getPyPlot_cMap('rainbow', [], [], '"/Users/hzhu/.pyenv/versions/3.11.13/bin/python"');
% colormap(cmp)

% exportgraphics(figure(1),'axissymslap_lap_dlpn0th_convergence.png','Resolution',200)
% exportgraphics(figure(2),'axissymslap_lap_dlpn0th_error.png','Resolution',200)

function s=mkcurve(Z,p,np)
s.p=p; s.Z=Z; s.tpan=linspace(0,pi,np+1)'; s=quadr(s,[],'p','G');
if ~isfield(s,'np'), s.np=numel(s.tpan)-1; end
end

function u = lapring(t,ysrc)
% axisymmetric (0th-mode) Laplace potential of ring charges via the elliptic-K ring formula:
%   Phi(rho,z) = sum_j q_j/(2 pi^2) K(k_j^2)/sqrt((rho+a_j)^2+(z-b_j)^2),  k_j^2 = 4 rho a_j/(...).
rho=real(t.x(:)); z=imag(t.x(:)); u=zeros(numel(rho),1);
for j=1:numel(ysrc.x)
  a=real(ysrc.x(j)); b=imag(ysrc.x(j));
  d2=(rho+a).^2+(z-b).^2; k2=4*rho*a./d2;
  u=u+ysrc.q(j)/(2*pi^2)*ellipke(k2)./sqrt(d2);
end
end

function g = lapringn(t,ysrc)
% n.grad of the axisymmetric ring potential, via azimuthal quadrature of the 3D point-charge gradient.
rho=real(t.x(:)); z=imag(t.x(:)); nr=real(t.nx(:)); nz=imag(t.nx(:));
Nphi=512; phi=2*pi*(0:Nphi-1)/Nphi; cphi=cos(phi); sphi=sin(phi);
g=zeros(numel(rho),1);
for j=1:numel(ysrc.x)
  a=real(ysrc.x(j)); b=imag(ysrc.x(j)); qj=ysrc.q(j);
  for ip=1:Nphi
    dx=rho-a*cphi(ip); dy=-a*sphi(ip); dz=z-b; r2=dx.^2+dy.^2+dz.^2;
    g=g - qj/Nphi*(nr.*dx + nz.*dz)./(4*pi*r2.^1.5);
  end
end
end
