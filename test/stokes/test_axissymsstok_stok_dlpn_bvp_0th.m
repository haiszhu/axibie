clearvars; close all; format short e;
addpath('/Users/hzhu/Documents/Github/axibie/utils');
addpath('/Users/hzhu/Documents/Github/axibie/alpertkernels');  % AxiKernel, AxiKernelT
addpath('/Users/hzhu/Documents/Github/axibie/matlab');

p=16; iside=1; iclosed=0; mu=1;
lam=0.75; Z=@(t) -(1.5+cos(t)).*(-sin(lam*pi*sin(t))+1i*cos(lam*pi*sin(t)));  % c-shape (open, poles t=0,pi)
Np=6:2:36;
ysrc.x=[0.70-1.20i; 0.50-1.80i; 1.30+0.20i]; ysrc.ws=ones(size(ysrc.x));     % ring forces INSIDE
F=[1;-0.7;0.5; 0.4;0.9;-0.6];                          % [F_rho(3); F_zeta(3)]

% 1. exact axi-Stokeslet: velocity (uexact) and traction sigma.n (texact, needs tt.nx)
uexact=@(tt) AxiStokeslet(ysrc,tt)*F;                  % 2N -> [u_rho; u_zeta]
texact=@(tt) AxiStokesletT(ysrc,tt)*F;                 % 2N -> [t_rho; t_zeta] (sigma.n)

% 2. exterior meridian target grid (rho>=0), exact velocity on it
nx=70; gx=(1:nx)/nx*4; ny=80; gy=((1:ny)/ny*2-1)*3;
[xx,yy]=meshgrid(gx,gy); zz=xx+1i*yy;
s0=mkcurve(Z,p,Np(end));
[IN,ON]=inpolygon(real(zz),imag(zz),real(s0.x),imag(s0.x)); ii=(~IN)&(~ON);
t=[]; t.x=zz(ii(:));
ue=uexact(t); uf=nan*(1+1i)*zz; uf(ii)=ue(1:end/2)+1i*ue(end/2+1:end);

% 3. h-refinement: combined-field Neumann solve (D'+S') tau = sigma.n (cond ~1e3, vs S'-only ~1e11),
%    eval velocity (D+S)[tau]
err=nan(size(Np));
for k=1:numel(Np)
  s=mkcurve(Z,p,Np(k)); N=numel(s.x); np=s.np;
  sx=s.x(:); snx=s.nx(:); sws=s.ws(:); swxp=s.wxp(:);
  tpan=s.tpan(:); sxlo=s.Z(s.tpan(1:end-1)); sxlo=sxlo(:); sxhi=s.Z(s.tpan(2:end)); sxhi=sxhi(:);
  fp=texact(s);                                        % Neumann data = exact traction (surface normal)
  A = axss_dlpn_blockmat_mex(N, sx, snx, p,np, sx,snx,sws,swxp,tpan,sxlo,sxhi, iside,iclosed,mu, []) ...
      + axss_slpn_blockmat_mex(N, sx, snx, p,np, sx,snx,sws,swxp,tpan,sxlo,sxhi, iside,iclosed,mu, []);
  tau=A\fp;
  Ae = axss_dlp_blockmat_mex(numel(t.x), t.x(:), p,np, sx,snx,sws,swxp,tpan,sxlo,sxhi, iside,iclosed,mu, []) ...
       + axss_slp_blockmat_mex(numel(t.x), t.x(:), p,np, sx,snx,sws,swxp,tpan,sxlo,sxhi, iside,iclosed,mu, []);
  tt=Ae*tau;
  u=nan*(1+1i)*zz; u(ii)=tt(1:end/2)+1i*tt(end/2+1:end);
  err(k)=max(abs(u(:)-uf(:)));
  fprintf('Stokes combined (D''+S'') 0th-mode Neumann BVP (D'' mex): np=%2d, N=%d: max err vs exact = %.3e\n',Np(k),N,err(k));
end

% 4. meridian error map + convergence
figure(1); clf; semilogy(Np,err,'o-k'); xlabel('n_p'); ylabel('max err'); grid on;
title('Stokes combined (D''+S'') 0th-mode Neumann BVP (D'' mex): h-refinement');
figure(2); clf; imagesc(gx,gy,log10(abs(u-uf))); colorbar; hold on;
fill(real([s.x;s.x(1)]),imag([s.x;s.x(1)]),'w'); plot(s.x,'-k'); plot(ysrc.x,'r.','MarkerSize',12);
clim([-15 -8]); colormap('jet'); axis equal tight;
xlabel('\rho'); ylabel('z'); title('Stokes combined (D''+S'') 0th-mode: log_{10} err vs exact axi-Stokeslet velocity');

% exportgraphics(figure(1),'axissymsstok_stok_dlpn0th_convergence.png','Resolution',200)
% exportgraphics(figure(2),'axissymsstok_stok_dlpn0th_error.png','Resolution',200)

function s=mkcurve(Z,p,np)
s.p=p; s.Z=Z; s.tpan=linspace(0,pi,np+1)'; s=quadr(s,[],'p','G');
if ~isfield(s,'np'), s.np=numel(s.tpan)-1; end
end
function G=AxiStokeslet(s,t)   % axibie exact axisymmetric Stokeslet velocity (source s, target t)
X=[real(t.x);imag(t.x)]; Y=[real(s.x);imag(s.x)]; M=length(X)/2; N=length(Y)/2;
x1=X(1:M); x2=X(1+M:2*M); y1=Y(1:N); y2=Y(1+N:2*N); G=zeros(2*M,2*N);
for j=1:M
  Ker=AxiKernel(y1,y2,x1(j),x2(j));
  G(j,1:N)=((Ker(:,2)+Ker(:,3)))'; G(j,N+1:end)=(Ker(:,4))';
  G(j+M,1:N)=(Ker(:,5))'; G(j+M,N+1:end)=((Ker(:,1)+Ker(:,6)))';
end
W=0.5*ones(size(s.x))/8/pi; G=G.*[repmat(W',2*M,1),repmat(W',2*M,1)]; G=2*G;
end
function G=AxiStokesletT(s,t)  % axibie exact axi-Stokeslet TRACTION sigma.n (target t with normal t.nx)
X=[real(s.x);imag(s.x)]; m=length(X)/2; x1=X(1:m); x2=X(1+m:2*m); mt=numel(t.x);
G1=zeros(mt,m); G2=G1; G3=G1; G4=G1;
for ne=1:mt
  wt=ones(size(s.x));
  Ker=AxiKernelT(x1,x2,real(t.x(ne)),imag(t.x(ne)),real(t.nx(ne)),imag(t.nx(ne)));
  G1(ne,:)=(wt.*(Ker(:,1)))'; G2(ne,:)=(wt.*Ker(:,2))'; G3(ne,:)=(wt.*Ker(:,3))'; G4(ne,:)=(wt.*(Ker(:,4)))';
end
G=[G1 G2;G3 G4]*(-3/4/pi);
end
