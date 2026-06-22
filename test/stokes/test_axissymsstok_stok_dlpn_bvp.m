clearvars; close all; format short e;
addpath('/Users/hzhu/Documents/Github/axibie/utils');
addpath('/Users/hzhu/Documents/Github/AxiStokes3D/utils');
addpath('/Users/hzhu/Documents/Github/AxiStokes3D/matlab');
% Combined-field (D' + S') exterior-Neumann Stokes BVP, all-modes, h-refinement: solve the per-mode 3x3
% traction systems (D'(s,s,m)+S'(s,s,m)) tau_m = t_m for the density, eval velocity (D+S)[tau].  D' is the
% hypersingular double-layer traction; adding S' makes A well-conditioned (cond ~1e3, 2nd-kind-type).
% Fully Fortran mex: solve via axss_{dlpn,slpn}_blockmat_nmode_mex, eval via axss_{dlp,slp}_blockmat_nmode_mex
% (each call returns all modes 0..pmodes at once -> index (:,:,m+1) in the per-mode loop).

y1=[0.30;0;0.20]; F1=[1;-0.7;0.5]; y2=[0.32*cos(pi/4);0.32*sin(pi/4);-0.25]; F2=[0.4;0.9;-0.6];
uex=@(X) stk(X,y1,F1)+stk(X,y2,F2);
p=16; np_vals=2:16; errmax=nan(1,numel(np_vals));

for kk=1:numel(np_vals)
  np=np_vals(kk); pmodes=2*np; nmodes=2*pmodes+1; nang=nmodes; mu=1;

  % 1. geometry (sphere) + Neumann data: 3D traction sigma.n on Gamma -> per-mode (azimuthal FFT)
  s=[]; s.p=p; s.Z=@(t) sin(t)-1i*cos(t); s.tpan=linspace(0,pi,np+1)'; s=quadr(s,[],'p','G');
  N=numel(s.x);
  s3d=axisym_to_3d_quadrature(real(s.x),imag(s.x),s.ws,nang);
  nr=real(s.nx).'; nz=imag(s.nx).';
  tr=zeros(nang,N); tp=tr; tz=tr;
  for i=1:nang
    th=2*pi*(i-1)/nang; X=s3d.x(:,(i-1)*N+(1:N)); N3=[nr*cos(th); nr*sin(th); nz];
    T=stktrac(X,N3,y1,F1)+stktrac(X,N3,y2,F2);
    tr(i,:)=T(1,:)*cos(th)+T(2,:)*sin(th); tp(i,:)=-T(1,:)*sin(th)+T(2,:)*cos(th); tz(i,:)=T(3,:);
  end
  tmr=fftshift(fft(tr,nmodes,1)/nang,1); tmp=fftshift(fft(tp,nmodes,1)/nang,1); tmz=fftshift(fft(tz,nmodes,1)/nang,1);

  % 2. combined-field (D'+S') Neumann solve  (D'+S') tau_m = t_m  (Fortran mex, all modes at once)
  ax=@(blk) blk(N,s.x(:),s.nx(:),p,np,s.x(:),s.nx(:),s.ws(:),s.wxp(:),s.tpan(:),s.xlo(:),s.xhi(:),pmodes,1,0,mu,[]);
  Dpn=ax(@axss_dlpn_blockmat_nmode_mex); Spn=ax(@axss_slpn_blockmat_nmode_mex);
  tau=cell(pmodes+1,1);
  for m=0:pmodes
    A=Dpn(:,:,m+1)+Spn(:,:,m+1);
    rhs=[tmr(m+pmodes+1,:).'; tmp(m+pmodes+1,:).'; tmz(m+pmodes+1,:).'];
    tau{m+1}=A\rhs;
  end

  % 3. 3d target eval: velocity (D+S)[tau]
  ng3=16; gv=linspace(-1.9,1.9,ng3); [GX,GY,GZ]=meshgrid(gv,gv,gv);
  rr=sqrt(GX(:).^2+GY(:).^2+GZ(:).^2); ext3=rr>1+1e-2;
  P=[GX(ext3).';GY(ext3).';GZ(ext3).']; M3=size(P,2);
  rho3=hypot(P(1,:),P(2,:)); z3=P(3,:); phi3=atan2(P(2,:),P(1,:));
  t3=[]; t3.x=(rho3+1i*z3).';
  ev=@(blk) blk(M3,t3.x(:),p,np,s.x(:),s.nx(:),s.ws(:),s.wxp(:),s.tpan(:),s.xlo(:),s.xhi(:),pmodes,1,0,mu,[]);
  Dev=ev(@axss_dlp_blockmat_nmode_mex); Sev=ev(@axss_slp_blockmat_nmode_mex);
  v3r=zeros(pmodes+1,M3); v3p=v3r; v3z=v3r;
  for m=0:pmodes
    v=(Dev(:,:,m+1)+Sev(:,:,m+1))*tau{m+1};
    v3r(m+1,:)=v(1:M3).'; v3p(m+1,:)=v(M3+1:2*M3).'; v3z(m+1,:)=v(2*M3+1:3*M3).';
  end
  ur3=real(v3r(1,:)); up3=real(v3p(1,:)); uz3=real(v3z(1,:));
  for m=1:pmodes
    e=exp(1i*m*phi3);
    ur3=ur3+2*real(e.*v3r(m+1,:)); up3=up3+2*real(e.*v3p(m+1,:)); uz3=uz3+2*real(e.*v3z(m+1,:));
  end
  ux3=ur3.*cos(phi3)-up3.*sin(phi3); uy3=ur3.*sin(phi3)+up3.*cos(phi3);

  % 4. log10 error
  Uex3=uex(P); err3=vecnorm([ux3;uy3;uz3]-Uex3)/max(vecnorm(Uex3));
  d3=rr(ext3).'-1; inb=d3<0.1; errmax(kk)=max(err3);
  fprintf('Stokes combined (D''+S'') all-modes Neumann BVP (Fortran mex): N=%d, np=%d, pmodes=%d, M3=%d\n',N,np,pmodes,M3);
  if any(inb), fprintf('  near band (d<0.1):  max err = %.3e  (%d pts)\n',max(err3(inb)),nnz(inb)); end
  fprintf('  far field (d>=0.1): max err = %.3e  (%d pts)\n',max(err3(~inb)),nnz(~inb));
end

figure(1),clf; semilogy(np_vals,errmax,'o-k'); xlabel('n_p'); ylabel('max err'); grid on;
title('Stokes combined (D''+S'') all-modes Neumann BVP: h-refinement, p_{modes}=2 n_p');

% scatter3 log10 error (last refinement)
figure(2),clf; hold on;
thm=linspace(0,2*pi,49);
Lx=real(s.x)*cos(thm); Ly=real(s.x)*sin(thm); Lz=imag(s.x)*ones(size(thm));
mesh(Lx,Ly,Lz,'EdgeColor',[0.55 0.55 0.55],'FaceColor','none','EdgeAlpha',0.35);
scatter3(P(1,:),P(2,:),P(3,:),18,log10(max(err3,1e-17)),'filled');
plot3(y1(1),y1(2),y1(3),'p','MarkerSize',16,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
plot3(y2(1),y2(2),y2(3),'p','MarkerSize',16,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
cb=colorbar; cb.Label.String='log_{10}|u_h-u_{exact}|/max|u_{exact}|';
clim([-12 -6]); axis equal; view(35,18); grid on; colormap('jet');
xlabel('x'); ylabel('y'); zlabel('z'); title('Stokes combined (D''+S'') all-modes Neumann: 3D target-grid velocity error');

% exportgraphics(figure(1),'axissymsstok_stok_dlpn_convergence.png','Resolution',200)
% exportgraphics(figure(2),'axissymsstok_stok_dlpn_error.png','Resolution',200)

function U=stk(X,y,F)
d=X-y; r=vecnorm(d); U=(1/(8*pi))*(F./r+(F.'*d).*d./r.^3);
end

function T=stktrac(X,n,y,F)
d=X-y; r=vecnorm(d); rF=F.'*d; rn=sum(n.*d,1); T=-(3/(4*pi))*(rF.*rn).*d./r.^5;
end
