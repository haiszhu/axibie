clearvars; format short e;
addpath('../../utils');
addpath('../../matlab');

% interior point charges -> exact exterior harmonic potential u_ex (scalar)
y1=[0.30;0;0.20]; q1=1.0; y2=[0.32*cos(pi/4);0.32*sin(pi/4);-0.25]; q2=-0.7;
uex=@(X) q1./(4*pi*vecnorm(X-y1)) + q2./(4*pi*vecnorm(X-y2));
p=16; np_vals=2:16; errmax=nan(1,numel(np_vals));

for kk=1:numel(np_vals)
  np=np_vals(kk); pmodes=2*np; nmodes=2*pmodes+1; nang=nmodes; iside=1; iclosed=0;

  % 1. geometry + Dirichlet data setup
  s=[]; s.p=p; s.Z=@(t) sin(t)-1i*cos(t); s.tpan=linspace(0,pi,np+1)'; s=quadr(s,[],'p','G');
  N=numel(s.x);
  sx=s.x(:); snx=s.nx(:); sws=s.ws(:); swxp=s.wxp(:);
  tpan=s.tpan(:); sxlo=s.Z(s.tpan(1:end-1)); sxlo=sxlo(:); sxhi=s.Z(s.tpan(2:end)); sxhi=sxhi(:);
  s3d=axisym_to_3d_quadrature(real(s.x),imag(s.x),s.ws,nang);
  f=zeros(nang,N);
  for i=1:nang
    f(i,:)=uex(s3d.x(:,(i-1)*N+(1:N)));
  end
  fm=fftshift(fft(f,nmodes,1)/nang,1);

  % 2. self modal matrix (all modes, scalar real block)
  Aself=axls_slp_blockmat_nmode_mex(N,sx,p,np,sx,snx,sws,swxp,tpan,sxlo,sxhi,pmodes,iside,iclosed,[]);

  % 3. solve  S_m[sigma_m] = u_m
  sig=cell(pmodes+1,1);
  for m=0:pmodes
    sig{m+1}=lsqminnorm(Aself(:,:,m+1), fm(m+pmodes+1,:).');
  end

  % 4. 3d target eval (single-layer potential S[sigma])
  ng3=16; gv=linspace(-1.9,1.9,ng3); [GX,GY,GZ]=meshgrid(gv,gv,gv);
  rr=sqrt(GX(:).^2+GY(:).^2+GZ(:).^2); ext3=rr>1+1e-2;
  P=[GX(ext3).';GY(ext3).';GZ(ext3).']; M3=size(P,2);
  rho3=hypot(P(1,:),P(2,:)); z3=P(3,:); phi3=atan2(P(2,:),P(1,:));
  Ae=axls_slp_blockmat_nmode_mex(M3,(rho3+1i*z3).',p,np,sx,snx,sws,swxp,tpan,sxlo,sxhi,pmodes,iside,iclosed,[]);
  v3=zeros(pmodes+1,M3);
  for m=0:pmodes
    v3(m+1,:)=(Ae(:,:,m+1)*sig{m+1}).';
  end
  u3=real(v3(1,:));
  for m=1:pmodes, u3=u3+2*real(exp(1i*m*phi3).*v3(m+1,:)); end

  % 5. log10 error
  Uex3=uex(P); err3=abs(u3-Uex3)/max(abs(Uex3));
  d3=rr(ext3).'-1; inb=d3<0.1; errmax(kk)=max(err3);
  fprintf('Laplace SLP all-modes Dirichlet BVP (delta mex): N=%d, np=%d, pmodes=%d, M3=%d\n',N,np,pmodes,M3);
  if any(inb), fprintf('  near band (d<0.1):  max err = %.3e  (%d pts)\n',max(err3(inb)),nnz(inb)); end
  fprintf('  far field (d>=0.1): max err = %.3e  (%d pts)\n',max(err3(~inb)),nnz(~inb));
end

figure(1),clf; semilogy(np_vals,errmax,'o-k'); xlabel('n_p'); ylabel('max err'); grid on;
title('Laplace SLP all-modes Dirichlet BVP (delta mex): h-refinement, p_{modes}=2 n_p');

% scatter3 log10 error (last refinement)
figure(2),clf; hold on;
thm=linspace(0,2*pi,49);
Lx=real(s.x)*cos(thm); Ly=real(s.x)*sin(thm); Lz=imag(s.x)*ones(size(thm));
mesh(Lx,Ly,Lz,'EdgeColor',[0.55 0.55 0.55],'FaceColor','none','EdgeAlpha',0.35);
scatter3(P(1,:),P(2,:),P(3,:),18,log10(max(err3,1e-17)),'filled');
plot3(y1(1),y1(2),y1(3),'p','MarkerSize',16,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
plot3(y2(1),y2(2),y2(3),'p','MarkerSize',16,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
cb=colorbar; cb.Label.String='log_{10}|u_h-u_{exact}|/max|u_{exact}|';
clim([-16 -12]); axis equal; view(35,18); grid on; colormap('jet');
xlabel('x'); ylabel('y'); zlabel('z'); title('Laplace SLP all-modes: 3D target-grid potential error');

% cmp = getPyPlot_cMap('rainbow', [], [], '"/Users/hzhu/.pyenv/versions/3.11.13/bin/python"');
% colormap(cmp)

% exportgraphics(figure(1),'axissymslap_lap_slp_convergence.png','Resolution',200)
% exportgraphics(figure(2),'axissymslap_lap_slp_error.png','Resolution',200)
