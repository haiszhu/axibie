function test_axissymstok_dirichlet
% Stokes dirichlet, MULTI-MODE -- combined-field (SLP+DLP) interior BVP.
%
clearvars; close all;
here=fileparts(mfilename('fullpath'));
addpath('/Users/hzhu/Documents/Github/chunkie/chunkie');     % +lege (lege.exps/pols, panel setup)
addpath('/Users/hzhu/Documents/Github/axibie/matlab');  % Fortran mex (axm_specialquad_*, axa_kernel_*)

% ---- geometry: torus tube = circle of radius R about (rc,0) ----
rc=3; R=1; mu=1.0;
Z =@(t)(rc+R*cos(t))+1i*(R*sin(t)); Zp=@(t)R*(-sin(t)+1i*cos(t));
Npan=16; p=16; [gx,gw]=lege.exps(p); edges=linspace(0,2*pi,Npan+1).';
ns=Npan*p; yk=zeros(ns,1); nvk=zeros(ns,1); wl=zeros(ns,1);
for kc=1:Npan
  a0=edges(kc); b0=edges(kc+1); hm=(b0-a0)/2; mid=(a0+b0)/2; idx=(kc-1)*p+(1:p); tt=mid+hm*gx;
  yp=Zp(tt); tau=yp./abs(yp); yk(idx)=Z(tt); nvk(idx)=-1i*tau; wl(idx)=abs(yp).*gw*hm;
end
rs=real(yk).'; zs=imag(yk).'; nrs=real(nvk).'; nzs=imag(nvk).';

src_pts=[ 1.0  -1.0   0.6;
          0.5   1.0  -0.8;
          2.5  -2.8   3.2];                 % source locations (columns)
src_frc=[ 1.0   2.0  -1.5;
          2.0  -1.0   0.7;
          3.0   1.0   2.0];                 % point forces
uex=@(X) stokeslet_field(X, src_pts, src_frc);

% ---- Dirichlet data f_m: sample azimuths, FFT to the modal meridional frame ----
Mtrunc=24; n_angles=2*Mtrunc+1; modes=-Mtrunc:Mtrunc; n_modes=numel(modes);
f3=zeros(n_angles,3*ns);
for j=1:n_angles
  th=(j-1)*2*pi/n_angles;
  X=[rs*cos(th); rs*sin(th); zs];
  f3(j,:)=cart_to_cyl_interleaved(uex(X), th);
end
f_m=fftshift(fft(f3,n_angles,1)/n_angles,1);

% ---- per-panel close-eval quadrature: be=2 upsampled (q=2p) + q->p Legendre projection ----
be=2; q=be*p; [gx2,gw2]=lege.exps(q); [~,~,uk]=lege.exps(p); Pmat=(lege.pols(gx2,p-1)).'*uk;
even_mask=[1 0 1;0 1 0;1 0 1]; odd_mask=[0 1 0;1 0 1;0 1 0];
zt_surf=(rs+1i*zs).';

% build the combined operator A_m = -I/2 + V^S_m + V^D_m, then solve A_m mu_m = f_m
A=cell(n_modes,1); for i=1:n_modes, A{i}=zeros(3*ns,3*ns); end
tic
for kc=1:Npan
  a0=edges(kc); b0=edges(kc+1); hm=(b0-a0)/2; mid=(a0+b0)/2; tt=mid+hm*gx;
  yp=Zp(tt); tau=yp./abs(yp); nv=-1i*tau; panlen=sum(abs(yp).*gw)*hm; pc=Z(mid); wlp=abs(yp).*gw*hm;
  cols=(kc-1)*p+(1:p); idx3=reshape(3*(cols-1)+(1:3).',[],1); w3p=reshape(repmat(wlp(:).',3,1),1,[]);
  sr=[real(Z(tt)).';imag(Z(tt)).']; sn=[real(nv).';imag(nv).'];
  yU=Z(mid+hm*gx2); ypU=Zp(mid+hm*gx2); za=Z(a0); zb=Z(b0);
  dpc=sqrt((rs-real(pc)).^2+(zs-imag(pc)).^2); nb=dpc<1.5*panlen; ni=find(nb); fi=find(~nb);
  if ~isempty(ni)
    nn=numel(ni); zn=zt_surf(ni); rows=reshape(3*(ni(:).'-1)+(1:3).',[],1);
    WS=reshape(axm_specialquad_slp_mex(nn,p,Mtrunc,mu,hm,0, real(zn),imag(zn), real(yU),imag(yU), real(ypU),imag(ypU), real(za),imag(za), real(zb),imag(zb), gw2, Pmat, []), Mtrunc+1,3*nn,3*p);
    WD=reshape(axm_specialquad_dlp_mex(nn,p,Mtrunc,mu,hm,0, real(zn),imag(zn), real(yU),imag(yU), real(ypU),imag(ypU), real(za),imag(za), real(zb),imag(zb), gw2, Pmat, []), Mtrunc+1,3*nn,3*p);
    for i=1:n_modes
      A{i}(rows,idx3)=squeeze(WS(abs(modes(i))+1,:,:))+squeeze(WD(abs(modes(i))+1,:,:));
    end
  end
  if ~isempty(fi)
    nf=numel(fi); rows=reshape(3*(fi(:).'-1)+(1:3).',[],1); tgr=rs(fi); tgz=zs(fi);
    for i=1:n_modes
      m=abs(modes(i));
      KS=axa_kernel_slp_mex(nf,p,sr(1,:),sr(2,:),tgr,tgz,m,mu,[]).*repmat(w3p,3*nf,1);
      KD=axa_kernel_dlp_mex(nf,p,sr(1,:),sr(2,:),sn(1,:),sn(2,:),tgr,tgz,m,mu,[]).*repmat(w3p,3*nf,1);
      A{i}(rows,idx3)=KS+KD;
    end
  end
end
% for i = 1:n_modes
%   m = abs(modes(i));
%   S = kernel();
%   S.eval = @(s,t) axa_kernel_slp_mex(size(t.r,2), size(s.r,2), ...
%                      s.r(1,:), s.r(2,:), t.r(1,:), t.r(2,:), m, mu, []);
%   S.opdims = [3,3]; S.sing = 'log';
%   D = kernel();
%   D.eval = @(s,t) axa_kernel_dlp_mex(size(t.r,2), size(s.r,2), ...
%                      s.r(1,:), s.r(2,:), t.r(1,:), t.r(2,:), m, mu, []);
%   D.opdims = [3,3]; D.sing = 'pv'; % is this true?
%   A{i} = chunkermat(chnkr, S, opts) + chunkermat(chnkr, D, opts);
% end
mu_m=zeros(n_modes,3*ns);
for i=1:n_modes
  sgn=sign(modes(i)); if sgn==0, sgn=1; end
  mask=[1,1i*sgn,1; 1i*sgn,1,1i*sgn; 1,1i*sgn,1];
  Ai=A{i}.*kron(ones(ns,ns),mask);
  mu_m(i,:)=gmres(Ai,f_m(i,:).',[],1e-12,3*ns).';
end
fprintf('\n test_axissymstok_dirichlet (combined SLP+DLP, Fortran close-eval mex, modes 0..%d)\n',Mtrunc);
fprintf('  operator build + solve: '); toc

% interior evaluation:  u = (V^S + V^D)[mu]  on a full 3D meshgrid inside the tube
xg=linspace(-4.5,4.5,31); yg=linspace(-4.5,4.5,31); zg=linspace(-1.2,1.2,11);
[Xg,Yg,Zg]=meshgrid(xg,yg,zg); xt=Xg(:).'; yt=Yg(:).'; zt3=Zg(:).';
rt=hypot(xt,yt); tht=atan2(yt,xt);
inside=((rt-rc).^2+zt3.^2)<R^2; gi=find(inside); Ng=numel(gi);
trg=rt(gi); tzg=zt3(gi); thg=tht(gi); ztg=(trg+1i*tzg).';   % meridian (rho,z) + azimuth per interior target

UE=cell(n_modes,1); UO=cell(n_modes,1);
for i=1:n_modes, UE{i}=zeros(3,Ng); UO{i}=zeros(3,Ng); end
for kc=1:Npan
  a0=edges(kc); b0=edges(kc+1); hm=(b0-a0)/2; mid=(a0+b0)/2; tt=mid+hm*gx;
  yp=Zp(tt); tau=yp./abs(yp); nv=-1i*tau; panlen=sum(abs(yp).*gw)*hm; pc=Z(mid); wlp=abs(yp).*gw*hm;
  cols=(kc-1)*p+(1:p); idx3=reshape(3*(cols-1)+(1:3).',[],1); w3p=reshape(repmat(wlp(:).',3,1),1,[]);
  sr=[real(Z(tt)).';imag(Z(tt)).']; sn=[real(nv).';imag(nv).'];
  yU=Z(mid+hm*gx2); ypU=Zp(mid+hm*gx2); za=Z(a0); zb=Z(b0);
  dpc=abs(ztg-pc).'; nb=dpc<1.5*panlen; ni=find(nb); fi=find(~nb);
  if ~isempty(ni)
    nn=numel(ni); zn=ztg(ni); EVn=kron(ones(nn,p),even_mask); ODn=kron(ones(nn,p),odd_mask);
    WS=reshape(axm_specialquad_slp_mex(nn,p,Mtrunc,mu,hm,0, real(zn),imag(zn), real(yU),imag(yU), real(ypU),imag(ypU), real(za),imag(za), real(zb),imag(zb), gw2, Pmat, []), Mtrunc+1,3*nn,3*p);
    WD=reshape(axm_specialquad_dlp_mex(nn,p,Mtrunc,mu,hm,0, real(zn),imag(zn), real(yU),imag(yU), real(ypU),imag(ypU), real(za),imag(za), real(zb),imag(zb), gw2, Pmat, []), Mtrunc+1,3*nn,3*p);
    for i=1:n_modes
      blk=squeeze(WS(abs(modes(i))+1,:,:))+squeeze(WD(abs(modes(i))+1,:,:)); dmu=mu_m(i,idx3).';
      UE{i}(:,ni)=UE{i}(:,ni)+reshape((blk.*EVn)*dmu,3,nn);
      UO{i}(:,ni)=UO{i}(:,ni)+reshape((blk.*ODn)*dmu,3,nn);
    end
  end
  if ~isempty(fi)
    nf=numel(fi); EVf=kron(ones(nf,p),even_mask); ODf=kron(ones(nf,p),odd_mask); tgr=trg(fi); tgz=tzg(fi);
    for i=1:n_modes
      m=abs(modes(i));
      KS=axa_kernel_slp_mex(nf,p,sr(1,:),sr(2,:),tgr,tgz,m,mu,[]).*repmat(w3p,3*nf,1);
      KD=axa_kernel_dlp_mex(nf,p,sr(1,:),sr(2,:),sn(1,:),sn(2,:),tgr,tgz,m,mu,[]).*repmat(w3p,3*nf,1);
      blk=KS+KD; dmu=mu_m(i,idx3).';
      UE{i}(:,fi)=UE{i}(:,fi)+reshape((blk.*EVf)*dmu,3,nf);
      UO{i}(:,fi)=UO{i}(:,fi)+reshape((blk.*ODf)*dmu,3,nf);
    end
  end
end
u_cyl=zeros(3,Ng);
for i=1:n_modes
  sgn=sign(modes(i)); if sgn==0, sgn=1; end
  u_cyl=u_cyl+real((UE{i}+1i*sgn*UO{i}).*exp(1i*modes(i)*thg));
end
% meridional (e_r,e_phi,e_z) -> Cartesian per target azimuth
u_sol=[u_cyl(1,:).*cos(thg)-u_cyl(2,:).*sin(thg);
       u_cyl(1,:).*sin(thg)+u_cyl(2,:).*cos(thg);
       u_cyl(3,:)];
u_true=uex([xt(gi); yt(gi); zt3(gi)]);
errg=vecnorm(u_sol-u_true,2,1);
dsurf=abs(sqrt((trg-rc).^2+tzg.^2)-R);
fprintf('  interior 3D grid: %d targets, max ||u-u_exact|| = %.3e\n', Ng, max(errg));
fprintf('    near-surface band (dist<0.05): max err = %.3e\n', max(errg(dsurf<0.05)));

% ---- figure: (1) Dirichlet boundary data on the torus surface, (2) interior log10 error ----
th_all=(0:n_angles-1)*2*pi/n_angles;
Sx=[reshape(rs(:)*cos(th_all),1,[]); reshape(rs(:)*sin(th_all),1,[]); reshape(repmat(zs(:),1,n_angles),1,[])];
u_surf=uex(Sx); spd=vecnorm(u_surf,2,1); sk=1:15:size(Sx,2);
figure('Position',[50 60 1340 660]);
ax1=subplot(1,2,1);
scatter3(Sx(1,:),Sx(2,:),Sx(3,:),16,spd,'filled'); hold on
quiver3(Sx(1,sk),Sx(2,sk),Sx(3,sk), u_surf(1,sk),u_surf(2,sk),u_surf(3,sk), 3,'k','LineWidth',0.4);
plot3(src_pts(1,:),src_pts(2,:),src_pts(3,:),'r^','MarkerFaceColor','r','MarkerSize',10);
axis equal; grid on; box on; view(35,20); colorbar; colormap(ax1,parula);
xlabel('x'); ylabel('y'); zlabel('z');
title('Dirichlet boundary data:  u on the torus surface   (color = |u|, arrows = u)');
legend({'surface nodes (|u|)','u (subsampled)','exterior Stokeslets'},'Location','northeast');
ax2=subplot(1,2,2);
scatter3(xt(gi),yt(gi),zt3(gi),26,log10(errg+1e-16),'filled');
axis equal; grid on; box on; view(35,20); colorbar; colormap(ax2,jet); caxis(ax2,[-16 -8]);
xlabel('x'); ylabel('y'); zlabel('z');
title(sprintf('log_{10} ||u_h - u_{exact}||   interior (combined SLP+DLP, modes 0..%d)',Mtrunc));
outpng=fullfile(here,'dirichlet_combined_fortran_mex.png');
saveas(gcf,outpng); fprintf('  saved -> %s\n',outpng);
end

% ===========================================================================
function u = stokeslet_field(X, pts, frc)
u=zeros(size(X));
for k=1:size(pts,2), u=u+stokeslet_velocity(X, pts(:,k), frc(:,k)); end
end

% ===========================================================================
function v = cart_to_cyl_interleaved(v_cart, theta)
n=size(v_cart,2); v=zeros(1,3*n);
v(1:3:end)= v_cart(1,:)*cos(theta)+v_cart(2,:)*sin(theta);
v(2:3:end)=-v_cart(1,:)*sin(theta)+v_cart(2,:)*cos(theta);
v(3:3:end)= v_cart(3,:);
end

% ===========================================================================
function u = stokeslet_velocity(x, y, f)
r=x-y; rmag=sqrt(sum(r.^2,1)); rdotf=f.'*r; u=(1/(8*pi)).*( f./rmag + r.*(rdotf./rmag.^3) );
end
