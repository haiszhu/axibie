clearvars; close all; format short e;
addpath('../../utils');
addpath('../../matlab');

y1=[0.30;0;0.20]; F1=[1;-0.7;0.5]; y2=[0.32*cos(pi/4);0.32*sin(pi/4);-0.25]; F2=[0.4;0.9;-0.6];
uex=@(X) stk(X,y1,F1)+stk(X,y2,F2);
p=16; np_vals=2:16; errmax=nan(1,numel(np_vals));

for kk=1:numel(np_vals)
  np=np_vals(kk); pmodes=2*np; nmodes=2*pmodes+1; nang=nmodes; iside=1; iclosed=0; mu=1;

  % 1. geometry + boundary data setup
  s=[]; s.p=p; s.Z=@(t) sin(t)-1i*cos(t); s.tpan=linspace(0,pi,np+1)'; s=quadr(s,[],'p','G');
  N=numel(s.x);
  sx=s.x(:); snx=s.nx(:); sws=s.ws(:); swxp=s.wxp(:);
  tpan=s.tpan(:); sxlo=s.Z(s.tpan(1:end-1)); sxlo=sxlo(:); sxhi=s.Z(s.tpan(2:end)); sxhi=sxhi(:);
  s3d=axisym_to_3d_quadrature(real(s.x),imag(s.x),s.ws,nang);
  fr=zeros(nang,N); fp=fr; fz=fr;
  for i=1:nang
    th=2*pi*(i-1)/nang; U=uex(s3d.x(:,(i-1)*N+(1:N)));
    fr(i,:)=U(1,:)*cos(th)+U(2,:)*sin(th); fp(i,:)=-U(1,:)*sin(th)+U(2,:)*cos(th); fz(i,:)=U(3,:);
  end
  fmr=fftshift(fft(fr,nmodes,1)/nang,1); fmp=fftshift(fft(fp,nmodes,1)/nang,1); fmz=fftshift(fft(fz,nmodes,1)/nang,1);

  % 2. self modal matrices (D_+ + S)
  Adlp=axss_dlp_blockmat_nmode_mex(N,sx,p,np,sx,snx,sws,swxp,tpan,sxlo,sxhi,pmodes,iside,iclosed,mu,[]);
  Aslp=axss_slp_blockmat_nmode_mex(N,sx,p,np,sx,snx,sws,swxp,tpan,sxlo,sxhi,pmodes,iside,iclosed,mu,[]);

  % 3. solve
  sig=cell(pmodes+1,1);
  for m=0:pmodes
    rhs=[fmr(m+pmodes+1,:).'; fmp(m+pmodes+1,:).'; fmz(m+pmodes+1,:).'];
    sig{m+1}=lsqminnorm(Adlp(:,:,m+1)+Aslp(:,:,m+1),rhs);
  end

  % 4. 3d target eval
  ng3=16; gv=linspace(-1.9,1.9,ng3); [GX,GY,GZ]=meshgrid(gv,gv,gv);
  rr=sqrt(GX(:).^2+GY(:).^2+GZ(:).^2); ext3=rr>1+1e-2;
  P=[GX(ext3).';GY(ext3).';GZ(ext3).']; M3=size(P,2);
  rho3=hypot(P(1,:),P(2,:)); z3=P(3,:); phi3=atan2(P(2,:),P(1,:));
  Aedlp=axss_dlp_blockmat_nmode_mex(M3,(rho3+1i*z3).',p,np,sx,snx,sws,swxp,tpan,sxlo,sxhi,pmodes,iside,iclosed,mu,[]);
  Aeslp=axss_slp_blockmat_nmode_mex(M3,(rho3+1i*z3).',p,np,sx,snx,sws,swxp,tpan,sxlo,sxhi,pmodes,iside,iclosed,mu,[]);
  v3r=zeros(pmodes+1,M3); v3p=v3r; v3z=v3r;
  for m=0:pmodes
    v=(Aedlp(:,:,m+1)+Aeslp(:,:,m+1))*sig{m+1};
    v3r(m+1,:)=v(1:M3).'; v3p(m+1,:)=v(M3+1:2*M3).'; v3z(m+1,:)=v(2*M3+1:3*M3).';
  end
  ur3=real(v3r(1,:)); up3=real(v3p(1,:)); uz3=real(v3z(1,:));
  for m=1:pmodes
    e=exp(1i*m*phi3);
    ur3=ur3+2*real(e.*v3r(m+1,:)); up3=up3+2*real(e.*v3p(m+1,:)); uz3=uz3+2*real(e.*v3z(m+1,:));
  end
  ux3=ur3.*cos(phi3)-up3.*sin(phi3); uy3=ur3.*sin(phi3)+up3.*cos(phi3);

  % 5. log10 error
  Uex3=uex(P); err3=vecnorm([ux3;uy3;uz3]-Uex3)/max(vecnorm(Uex3));
  d3=rr(ext3).'-1; inb=d3<0.1; errmax(kk)=max(err3);
  fprintf('Stokes combined (D_+ + S) all-modes BVP (delta mex): N=%d, np=%d, pmodes=%d, M3=%d\n',N,np,pmodes,M3);
  if any(inb), fprintf('  near band (d<0.1):  max err = %.3e  (%d pts)\n',max(err3(inb)),nnz(inb)); end
  fprintf('  far field (d>=0.1): max err = %.3e  (%d pts)\n',max(err3(~inb)),nnz(~inb));
end

figure(1),clf; semilogy(np_vals,errmax,'o-k'); xlabel('n_p'); ylabel('max err'); grid on;
title('Stokes combined (D_+ + S) all-modes BVP: h-refinement, p_{modes}=2 n_p');

% scatter3 log10 error (last refinement)
figure(2),clf; hold on;
thm=linspace(0,2*pi,49);
Lx=real(s.x)*cos(thm); Ly=real(s.x)*sin(thm); Lz=imag(s.x)*ones(size(thm));
mesh(Lx,Ly,Lz,'EdgeColor',[0.55 0.55 0.55],'FaceColor','none','EdgeAlpha',0.35);
scatter3(P(1,:),P(2,:),P(3,:),18,log10(max(err3,1e-17)),'filled');
plot3(y1(1),y1(2),y1(3),'p','MarkerSize',16,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
plot3(y2(1),y2(2),y2(3),'p','MarkerSize',16,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
cb=colorbar; cb.Label.String='log_{10}|u_h-u_{exact}|/max|u_{exact}|';
clim([-15 -8]); axis equal; view(35,18); grid on; colormap('jet');
xlabel('x'); ylabel('y'); zlabel('z'); title('Stokes combined (D_+ + S) all-modes: 3D target-grid velocity error');

% exportgraphics(figure(1),'axissymsstok_stok_dlp_convergence.png','Resolution',200)
% exportgraphics(figure(2),'axissymsstok_stok_dlp_error.png','Resolution',200)

function U=stk(X,y,F)
d=X-y; r=vecnorm(d); U=(1/(8*pi))*(F./r+(F.'*d).*d./r.^3);
end
