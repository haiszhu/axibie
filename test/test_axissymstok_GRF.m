function test_axissymstok_GRF
% Stokes Green's representation, MULTI-MODE

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
rs=real(yk).'; zs=imag(yk).'; nrs=real(nvk).'; nzs=imag(nvk).'; src.r=[rs;zs]; src.n=[nrs;nzs];

% ---- manufactured multi-mode field: one OFF-AXIS exterior Stokeslet ----
pt_force=[5.0; 0.0; 0.5]; y_force=[1.0; 0.5; 0.7];     % off-axis (rho=5) -> all azimuthal modes

% ---- modal surface densities: sample n_angles azimuths, FFT, keep |m|<=Mtrunc ----
Mtrunc=8; n_angles=64; modes=(-floor(n_angles/2):ceil(n_angles/2)-1);
sig3=zeros(n_angles,3*ns); tau3=sig3;
for j=1:n_angles
  th=(j-1)*2*pi/n_angles;
  X =[rs*cos(th); rs*sin(th); zs];  Nx=[nrs*cos(th); nrs*sin(th); nzs];
  uc=stokeslet_velocity(X,y_force,pt_force); tc=stokeslet_traction(X,Nx,y_force,pt_force);
  sig3(j,:)=cart_to_cyl_interleaved(tc,th); tau3(j,:)=cart_to_cyl_interleaved(-uc,th);
end
sig_m=fftshift(fft(sig3,n_angles,1)/n_angles,1); tau_m=fftshift(fft(tau3,n_angles,1)/n_angles,1);

% ---- meridian meshgrid targets (theta=0) ----
nr=56; nz=56; rg=linspace(rc-2*R,rc+2*R,nr); zg=linspace(-2*R,2*R,nz);
[RR,ZZ]=meshgrid(rg,zg); ntg=numel(RR); tr=RR(:).'; tz=ZZ(:).'; tg.r=[tr;tz];
inside=((tr-rc).^2+tz.^2)<R^2;
u_ex=stokeslet_velocity([tr;zeros(1,ntg);tz],y_force,pt_force); u_exp=u_ex.*inside;  % exact (theta=0)

% ---- SUBSECTION 1: NAIVE multi-mode field, sum modes 0..Mtrunc (e^{i m 0}=1) ----
even_mask=[1 0 1;0 1 0;1 0 1]; odd_mask=[0 1 0;1 0 1;0 1 0];
EVEN=kron(ones(ntg,ns),even_mask); ODD=kron(ones(ntg,ns),odd_mask);
w3=reshape(repmat(wl(:).',3,1),1,[]); urep_naive=zeros(3*ntg,1);
for m=0:Mtrunc
  Sm=axa_kernel_slp_mex(ntg,ns,src.r(1,:),src.r(2,:),tg.r(1,:),tg.r(2,:),m,mu,[]).*repmat(w3,3*ntg,1);
  Dm=axa_kernel_dlp_mex(ntg,ns,src.r(1,:),src.r(2,:),src.n(1,:),src.n(2,:),tg.r(1,:),tg.r(2,:),m,mu,[]).*repmat(w3,3*ntg,1);
  ip=find(modes==m,1);
  if m==0
    urep_naive=urep_naive+(Sm.*EVEN)*sig_m(ip,:).'+(Dm.*EVEN)*tau_m(ip,:).';
  else
    im=find(modes==-m,1);
    sPp=sig_m(ip,:).'; sPm=sig_m(im,:).'; tPp=tau_m(ip,:).'; tPm=tau_m(im,:).';
    SE=Sm.*EVEN; SO=Sm.*ODD; DE=Dm.*EVEN; DO=Dm.*ODD;
    urep_naive=urep_naive+SE*(sPp+sPm)+1i*SO*(sPp-sPm)+DE*(tPp+tPm)+1i*DO*(tPp-tPm);
  end
end
u_rep=reshape(real(urep_naive),3,ntg); err=vecnorm(u_rep-u_exp,2,1);
dsurf=abs(sqrt((tr-rc).^2+tz.^2)-R); far=dsurf>0.3; near=dsurf<0.05; ext=~inside;
fprintf('\n test_axissymstok_GRF (Fortran close-eval mex, modes 0..%d, EXTERIOR targets only)\n',Mtrunc);
fprintf('  SUBSECTION 1 (naive):       exterior far  (>0.3)  max err = %.3e\n',max(err(far&ext)));
fprintf('  SUBSECTION 1 (naive):       exterior near (<0.05) max err = %.3e   (naive quadrature fails)\n',max(err(near&ext)));

% ---- SUBSECTION 2: near-target correction via the FORTRAN close-eval mex, all modes 0..Mtrunc ----
clear Sm Dm SE SO DE DO EVEN ODD
urep_cor=urep_naive; emask=even_mask; omask=odd_mask;
% panel-independent close-eval quadrature: be=2 upsampled weights + q->p Legendre projection
be=2; q=be*p; [gx2,gw2]=lege.exps(q); [~,~,uk]=lege.exps(p); Pmat=(lege.pols(gx2,p-1)).'*uk;
for kc=1:Npan
  a0=edges(kc); b0=edges(kc+1); hm=(b0-a0)/2; mid=(a0+b0)/2; tt=mid+hm*gx;
  yp=Zp(tt); tau=yp./abs(yp); nv=-1i*tau; panlen=sum(abs(yp).*gw)*hm; pc=Z(mid); wlp=abs(yp).*gw*hm;
  srcp.r=[real(Z(tt)).'; imag(Z(tt)).']; srcp.n=[real(nv).'; imag(nv).'];
  tt2=mid+hm*gx2; yU=Z(tt2); ypU=Zp(tt2); za=Z(a0); zb=Z(b0);   % upsampled (q=2p) panel for the operator mex
  cols=(kc-1)*p+(1:p); w3p=reshape(repmat(wlp(:).',3,1),1,[]);
  dpc=sqrt((tr-real(pc)).^2+(tz-imag(pc)).^2); nearp=dpc<1.5*panlen;
  for sd=['e','i']
    if sd=='e', sel=nearp & ~inside; else, sel=nearp & inside; end
    gi=find(sel); if isempty(gi), continue; end
    nn=numel(gi); zn=(tr(gi)+1i*tz(gi)).'; tgn.r=[tr(gi);tz(gi)];
    rows=reshape(3*(gi(:).'-1)+(1:3).',[],1); iside=double(sd=='e');
    EVENb=kron(ones(nn,p),emask); ODDb=kron(ones(nn,p),omask);
    WS=reshape(axm_specialquad_slp_mex(nn,p,Mtrunc,mu,hm,iside, real(zn),imag(zn), real(yU),imag(yU), real(ypU),imag(ypU), real(za),imag(za), real(zb),imag(zb), gw2, Pmat, []), Mtrunc+1, 3*nn, 3*p);
    WD=reshape(axm_specialquad_dlp_mex(nn,p,Mtrunc,mu,hm,iside, real(zn),imag(zn), real(yU),imag(yU), real(ypU),imag(ypU), real(za),imag(za), real(zb),imag(zb), gw2, Pmat, []), Mtrunc+1, 3*nn, 3*p);
    for m=0:Mtrunc
      KS=axa_kernel_slp_mex(nn,p,srcp.r(1,:),srcp.r(2,:),tgn.r(1,:),tgn.r(2,:),m,mu,[]).*repmat(w3p,3*nn,1);
      KD=axa_kernel_dlp_mex(nn,p,srcp.r(1,:),srcp.r(2,:),srcp.n(1,:),srcp.n(2,:),tgn.r(1,:),tgn.r(2,:),m,mu,[]).*repmat(w3p,3*nn,1);
      dSraw=reshape(WS(m+1,:,:),3*nn,3*p)-KS; dDraw=reshape(WD(m+1,:,:),3*nn,3*p)-KD;
      ip=find(modes==m,1); idx3=reshape(3*(cols-1)+(1:3).',[],1);
      if m==0
        urep_cor(rows)=urep_cor(rows)+(dSraw.*EVENb)*sig_m(ip,idx3).'+(dDraw.*EVENb)*tau_m(ip,idx3).';
      else
        im=find(modes==-m,1);
        sPp=sig_m(ip,idx3).'; sPm=sig_m(im,idx3).'; tPp=tau_m(ip,idx3).'; tPm=tau_m(im,idx3).';
        SEb=dSraw.*EVENb; SOb=dSraw.*ODDb; DEb=dDraw.*EVENb; DOb=dDraw.*ODDb;
        urep_cor(rows)=urep_cor(rows)+SEb*(sPp+sPm)+1i*SOb*(sPp-sPm)+DEb*(tPp+tPm)+1i*DOb*(tPp-tPm);
      end
    end
  end
end
ucor=reshape(real(urep_cor),3,ntg); err2=vecnorm(ucor-u_exp,2,1); band=dsurf<0.1;
fprintf('  SUBSECTION 2 (Fortran mex): exterior near (<0.1)  max err = %.3e   (close-eval floor, modes 0..8)\n',max(err2(band&ext)));
fprintf('  SUBSECTION 2 (Fortran mex): exterior far  (>0.3)  max err = %.3e\n',max(err2(far&ext)));

% ---- figures (EXTERIOR targets only; interior blanked) ----
extg=reshape(ext,nz,nr); Lna=reshape(log10(err+1e-16),nz,nr); Lco=reshape(log10(err2+1e-16),nz,nr);
figure('Position',[80 80 1280 540]); th=linspace(0,2*pi,200); cc=rc+R*cos(th); ss=R*sin(th);
subplot(1,2,1); h1=imagesc(rg,zg,Lna); set(h1,'AlphaData',extg); set(gca,'YDir','normal');
axis equal tight; colorbar; caxis([-16 0]); hold on; plot(cc,ss,'k-','LineWidth',1.5);
xlabel('\rho'); ylabel('z'); title('NAIVE kernels (subsection 1)');
subplot(1,2,2); h2=imagesc(rg,zg,Lco); set(h2,'AlphaData',extg); set(gca,'YDir','normal');
axis equal tight; colorbar; caxis([-16 0]); hold on; plot(cc,ss,'k-','LineWidth',1.5);
xlabel('\rho'); ylabel('z'); title('FORTRAN-mex kernel-split close-eval (subsection 2)');
sgtitle('log_{10} | u_h - 0 |  EXTERIOR targets  (Stokes Green repn, modes 0..8, Fortran close-eval mex)');
outpng=fullfile(here,'green_GRF_fortran_mex.png');
saveas(gcf,outpng); fprintf('  saved -> %s\n',outpng);
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

% ===========================================================================
function trac = stokeslet_traction(x, n, y, f)
r=x-y; rmag=sqrt(sum(r.^2,1)); rdotf=f.'*r; rdotn=sum(r.*n,1);
trac=-(3/(4*pi)).*r.*(rdotf.*rdotn./rmag.^5);
end
