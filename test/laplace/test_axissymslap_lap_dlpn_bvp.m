clearvars; format short e;
addpath('../../utils');
addpath('../../matlab');

% ---- geometry switch (mirror axibie/testStokesDir.m / the SLP, SLPn, DLP tests) ----
shape = 'cshape';                                   % 'sphere' | 'ellipse' | 'cshape'
p=16; iside=1; iclosed=0;
if strcmp(shape,'sphere')
  Z=@(t) sin(t)-1i*cos(t);                          % unit sphere
  np_vals=2:16;
  y1=[0.30;0;0.20]; q1=1.0; y2=[0.32*cos(pi/4);0.32*sin(pi/4);-0.25]; q2=-0.7;
  gv=linspace(-1.9,1.9,16);
elseif strcmp(shape,'ellipse')
  a=2; Z=@(t) a^(-1/3)*sin(t)-1i*a^(2/3)*cos(t);    % volume-conserving prolate ellipsoid
  np_vals=2:16;
  y1=[0.20;0;0.40]; q1=1.0; y2=[0.25*cos(pi/4);0.25*sin(pi/4);-0.50]; q2=-0.7;
  gv=linspace(-2.4,2.4,18);
elseif strcmp(shape,'cshape')
  lam=0.75; Z=@(t) -(1.5+cos(t)).*(-sin(lam*pi*sin(t))+1i*cos(lam*pi*sin(t)));  % c-shape (near-axis poles)
  np_vals=4:2:24;
  y1=[0.70;0;-1.20]; q1=1.0; y2=[1.30*cos(pi/4);1.30*sin(pi/4);0.20]; q2=-0.7;
  gv=linspace(-3,3,24);
end

% interior point charges -> exact exterior harmonic potential u_ex; Neumann data g = dn u_ex.
% Combined-field Neumann: represent u = (S + D)[sigma], solve (-1/2 I + S' + D') sigma = dn u_ex
% (D' = d_n d_n' G is hypersingular and continuous; it regularizes the single-layer Neumann resonances).
uex =@(X) q1./(4*pi*vecnorm(X-y1)) + q2./(4*pi*vecnorm(X-y2));
gradu=@(X) -q1*(X-y1)./(4*pi*vecnorm(X-y1).^3) - q2*(X-y2)./(4*pi*vecnorm(X-y2).^3);

% ---- exterior 3D target grid (once): outside the revolved meridian, >1e-2 off the surface ----
mer=Z(linspace(0,pi,400).');                        % fine meridian for inpolygon + distance
[GX,GY,GZ]=meshgrid(gv,gv,gv); rho3=hypot(GX(:),GY(:)); z3=GZ(:); zt=rho3+1i*z3;
dmer=inf(size(zt)); for j=1:numel(mer), dmer=min(dmer,abs(zt-mer(j))); end
ext3=(~inpolygon(rho3,z3,real(mer),imag(mer))) & (dmer>1e-2);
P=[GX(ext3).';GY(ext3).';GZ(ext3).']; M3=size(P,2);
rho3=hypot(P(1,:),P(2,:)); z3=P(3,:); phi3=atan2(P(2,:),P(1,:)); d3=dmer(ext3).';

errmax=nan(1,numel(np_vals));
for kk=1:numel(np_vals)
  np=np_vals(kk); pmodes=2*np; nmodes=2*pmodes+1; nang=nmodes;

  % 1. geometry + Neumann (dn u_ex) data
  s=[]; s.p=p; s.Z=Z; s.tpan=linspace(0,pi,np+1)'; s=quadr(s,[],'p','G');
  N=numel(s.x);
  sx=s.x(:); snx=s.nx(:); sws=s.ws(:); swxp=s.wxp(:); tpan=s.tpan(:); sxlo=s.xlo(:); sxhi=s.xhi(:);
  nr=real(s.nx).'; nz=imag(s.nx).'; g=zeros(nang,N);
  for i=1:nang
    th=2*pi*(i-1)/nang; Y=[real(s.x).'*cos(th); real(s.x).'*sin(th); imag(s.x).']; N3=[nr*cos(th);nr*sin(th);nz];
    g(i,:)=sum(N3.*gradu(Y),1);
  end
  gm=fftshift(fft(g,nmodes,1)/nang,1);

  % 2. self modal matrices (-1/2 I + S' + D')  (target normal = surface normal)
  % FROZEN worker path (uncomment to compare): Sp=axls_slpn_blockmat_nmode_mex(N,sx,snx,p,np,sx,snx,sws,swxp,tpan,sxlo,sxhi,pmodes,iside,iclosed,[]);
  Sp=real(axp_modemat_setup_mex(1,2,0,1,1,p,np,pmodes,iside,iclosed,[1;N+1],[1;np+2],N,np+1, ...
      sx,snx,sws,swxp,tpan,0,[1;1],zeros(3,0),zeros(3,0),eye(3),zeros(3,1),N,N));
  % FROZEN worker path (uncomment to compare): Dp=axls_dlpn_blockmat_nmode_mex(N,sx,snx,p,np,sx,snx,sws,swxp,tpan,sxlo,sxhi,pmodes,iside,iclosed,[]);
  Dp=real(axp_modemat_setup_mex(1,4,0,1,1,p,np,pmodes,iside,iclosed,[1;N+1],[1;np+2],N,np+1, ...
      sx,snx,sws,swxp,tpan,0,[1;1],zeros(3,0),zeros(3,0),eye(3),zeros(3,1),N,N));

  % 3. solve  (-1/2 I + S'_m + D'_m) sigma_m = dn u_m
  sig=cell(pmodes+1,1);
  for m=0:pmodes
    sig{m+1}=(Sp(:,:,m+1)+Dp(:,:,m+1))\gm(m+pmodes+1,:).';
  end

  % 4. 3d target eval (combined potential (S+D)[sigma])
  % FROZEN worker path (uncomment to compare): Se=axls_slp_blockmat_nmode_mex(M3,(rho3+1i*z3).',p,np,sx,snx,sws,swxp,tpan,sxlo,sxhi,pmodes,iside,iclosed,[]);
  Se=real(axp_modemat_setup_mex(1,1,0,3,1,p,np,pmodes,iside,iclosed,[1;N+1],[1;np+2],N,np+1, ...
      sx,snx,sws,swxp,tpan,M3,[1;M3+1],P,zeros(3,M3),eye(3),zeros(3,1),M3,N));
  % FROZEN worker path (uncomment to compare): De=axls_dlp_blockmat_nmode_mex(M3,(rho3+1i*z3).',p,np,sx,snx,sws,swxp,tpan,sxlo,sxhi,pmodes,iside,iclosed,[]);
  De=real(axp_modemat_setup_mex(1,3,0,3,1,p,np,pmodes,iside,iclosed,[1;N+1],[1;np+2],N,np+1, ...
      sx,snx,sws,swxp,tpan,M3,[1;M3+1],P,zeros(3,M3),eye(3),zeros(3,1),M3,N));
  v3=zeros(pmodes+1,M3);
  for m=0:pmodes, v3(m+1,:)=((Se(:,:,m+1)+De(:,:,m+1))*sig{m+1}).'; end
  u3=real(v3(1,:));
  for m=1:pmodes, u3=u3+2*real(exp(1i*m*phi3).*v3(m+1,:)); end

  % 5. log10 error
  Uex3=uex(P); err3=abs(u3-Uex3)/max(abs(Uex3));
  inb=d3<0.1; errmax(kk)=max(err3);
  fprintf('Laplace combined (S''+D'') all-modes Neumann BVP [%s] (delta mex): N=%d, np=%d, pmodes=%d, M3=%d\n',shape,N,np,pmodes,M3);
  if any(inb), fprintf('  near band (d<0.1):  max err = %.3e  (%d pts)\n',max(err3(inb)),nnz(inb)); end
  fprintf('  far field (d>=0.1): max err = %.3e  (%d pts)\n',max(err3(~inb)),nnz(~inb));
end

figure(1),clf; semilogy(np_vals,errmax,'o-k'); xlabel('n_p'); ylabel('max err'); grid on;
title(sprintf('Laplace combined (S''+D'') all-modes Neumann BVP [%s]: h-refinement, p_{modes}=2 n_p',shape));

% scatter3 log10 error (last refinement)
figure(2),clf; hold on;
thm=linspace(0,2*pi,49);
Lx=real(s.x)*cos(thm); Ly=real(s.x)*sin(thm); Lz=imag(s.x)*ones(size(thm));
mesh(Lx,Ly,Lz,'EdgeColor',[0.55 0.55 0.55],'FaceColor','none','EdgeAlpha',0.35);
scatter3(P(1,:),P(2,:),P(3,:),18,log10(max(err3,1e-17)),'filled');
plot3(y1(1),y1(2),y1(3),'p','MarkerSize',16,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
plot3(y2(1),y2(2),y2(3),'p','MarkerSize',16,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
cb=colorbar; cb.Label.String='log_{10}|u_h-u_{exact}|/max|u_{exact}|';
clim([-16 -11]); axis equal; view(35,18); grid on; colormap('jet');
xlabel('x'); ylabel('y'); zlabel('z'); title(sprintf('Laplace combined (S''+D'') all-modes Neumann [%s]: 3D target-grid potential error',shape));

% cmp = getPyPlot_cMap('rainbow', [], [], '"/Users/hzhu/.pyenv/versions/3.11.13/bin/python"');
% colormap(cmp)

% exportgraphics(figure(1),'axissymslap_lap_dlpn_convergence.png','Resolution',200)
% exportgraphics(figure(2),'axissymslap_lap_dlpn_error.png','Resolution',200)
