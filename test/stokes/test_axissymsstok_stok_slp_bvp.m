clearvars; format short e;
addpath('../../utils');
addpath('../../matlab');

% ---- geometry switch (mirror the laplace tests: sphere | ellipse | cshape) ----
shape = 'cshape';                                   % 'sphere' | 'ellipse' | 'cshape'
p=16;
if strcmp(shape,'sphere')
  Z=@(t) sin(t)-1i*cos(t); np_vals=2:16;
  y1=[0.30;0;0.20]; y2=[0.32*cos(pi/4);0.32*sin(pi/4);-0.25]; gv=linspace(-1.9,1.9,16);
elseif strcmp(shape,'ellipse')
  a=2; Z=@(t) a^(-1/3)*sin(t)-1i*a^(2/3)*cos(t); np_vals=2:16;
  y1=[0.20;0;0.40]; y2=[0.25*cos(pi/4);0.25*sin(pi/4);-0.50]; gv=linspace(-2.4,2.4,18);
elseif strcmp(shape,'cshape')
  lam=0.75; Z=@(t) -(1.5+cos(t)).*(-sin(lam*pi*sin(t))+1i*cos(lam*pi*sin(t))); np_vals=6:2:24;
  y1=[0.70;0;-1.20]; y2=[1.30*cos(pi/4);1.30*sin(pi/4);0.20]; gv=linspace(-3,3,24);
end
F1=[1;-0.7;0.5]; F2=[0.4;0.9;-0.6];                 % interior Stokeslet forces (shape-independent)
uex=@(X) stk(X,y1,F1)+stk(X,y2,F2);
errmax=nan(1,numel(np_vals));
% ---- exterior 3D target grid (once): outside the revolved meridian, >1e-2 off the surface ----
mer=Z(linspace(0,pi,400).');
[GX,GY,GZ]=meshgrid(gv,gv,gv); rhog=hypot(GX(:),GY(:)); zg=GZ(:); ztg=rhog+1i*zg;
dmer=inf(size(ztg)); for jmer=1:numel(mer), dmer=min(dmer,abs(ztg-mer(jmer))); end
ext3=(~inpolygon(rhog,zg,real(mer),imag(mer))) & (dmer>1e-2);
P=[GX(ext3).';GY(ext3).';GZ(ext3).']; M3=size(P,2);
rho3=hypot(P(1,:),P(2,:)); z3=P(3,:); phi3=atan2(P(2,:),P(1,:)); d3=dmer(ext3).';

% ---- kdtree -----
backend = 1;
leafsize = 64;
pts = [rho3(:), z3(:), zeros(numel(z3(:)),1)];
axk_kdtree_build_mex(backend, pts, leafsize);

for kk=1:numel(np_vals)
  np=np_vals(kk); pmodes=2*np; nmodes=2*pmodes+1; nang=nmodes; iside=1; iclosed=0; mu=1;

  % 1. geometry + boundary data setup
  s=[]; s.p=p; s.Z=Z; s.tpan=linspace(0,pi,np+1)'; s=quadr(s,[],'p','G');
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

  % 2. self modal matrix
  Aself=axss_slp_blockmat_nmode_mex(N,sx,p,np,sx,snx,sws,swxp,tpan,sxlo,sxhi,pmodes,iside,iclosed,mu,[]);

  % 3. solve
  sig=cell(pmodes+1,1);
  for m=0:pmodes
    rhs=[fmr(m+pmodes+1,:).'; fmp(m+pmodes+1,:).'; fmz(m+pmodes+1,:).'];
    sig{m+1}=lsqminnorm(Aself(:,:,m+1),rhs);
  end

  % 4. near target
  isnear = false(M3,1);
  for k = 1:np
    jj  = (k-1)*p + (1:p);
    pk  = s.x(jj);
    qr  = mean(real(pk));  qz = mean(imag(pk));
    qradii = 1.75*sum(s.ws(jj));
    idxcin = axk_kdtree_ball_mex([qr, qz, 0], qradii, M3);
    isnear(idxcin) = true;
  end

  % 5. 3d target eval
  % targets P/rho3/z3/phi3/M3/d3 are precomputed once above (shape-dependent grid)
  zt = (rho3+1i*z3).';
  Mn = nnz(isnear);   Mf = nnz(~isnear);
  ux3=zeros(1,M3); uy3=zeros(1,M3); uz3=zeros(1,M3);   % near + far 分别填进来
  % near
  Aen = axss_slp_blockmat_nmode_mex(Mn, zt(isnear),  p,np,sx,snx,sws,swxp,tpan,sxlo,sxhi,pmodes,iside,iclosed,mu,[]);
  vr=zeros(pmodes+1,Mn); vp=vr; vz=vr;
  for m=0:pmodes
    v=Aen(:,:,m+1)*sig{m+1};
    vr(m+1,:)=v(1:Mn).'; vp(m+1,:)=v(Mn+1:2*Mn).'; vz(m+1,:)=v(2*Mn+1:3*Mn).';
  end
  phin = phi3(isnear);                       % 1×Mn
  urn=real(vr(1,:)); upn=real(vp(1,:)); uzn=real(vz(1,:));
  for m=1:pmodes
    e=exp(1i*m*phin);
    urn=urn+2*real(e.*vr(m+1,:)); upn=upn+2*real(e.*vp(m+1,:)); uzn=uzn+2*real(e.*vz(m+1,:));
  end
  ux3(isnear)=urn.*cos(phin)-upn.*sin(phin);
  uy3(isnear)=urn.*sin(phin)+upn.*cos(phin);
  uz3(isnear)=uzn;
  % far
  SR = cell2mat(cellfun(@(c) c(1:N).',      sig,'uni',0));
  SP = cell2mat(cellfun(@(c) c(N+1:2*N).',  sig,'uni',0));
  SZ = cell2mat(cellfun(@(c) c(2*N+1:3*N).',sig,'uni',0));
  Fx=zeros(N,nang); Fy=Fx; Fz=Fx;
  for i=1:nang
    th=2*pi*(i-1)/nang;
    tr=real(SR(1,:)); tp=real(SP(1,:)); tz=real(SZ(1,:));
    for m=1:pmodes
      e=exp(1i*m*th);
      tr=tr+2*real(e*SR(m+1,:)); tp=tp+2*real(e*SP(m+1,:)); tz=tz+2*real(e*SZ(m+1,:));
    end
    Fx(:,i)=(tr*cos(th)-tp*sin(th)).';
    Fy(:,i)=(tr*sin(th)+tp*cos(th)).';
    Fz(:,i)=tz.';
  end
  % Afar = Sto3dSLPmat(struct('x',P(:,~isnear)), s3d);
  % ufar = Afar*[Fx(:); Fy(:); Fz(:)];
  % ux3(~isnear)=ufar(1:Mf).'; uy3(~isnear)=ufar(Mf+1:2*Mf).'; uz3(~isnear)=ufar(2*Mf+1:3*Mf).';
  ufar = axt_sto3dslp_eval_mex(P(:,~isnear), s3d.x, s3d.w, [Fx(:),Fy(:),Fz(:)]);
  ux3(~isnear)=ufar(1,:); uy3(~isnear)=ufar(2,:); uz3(~isnear)=ufar(3,:);

  % 6. log10 error
  Uex3=uex(P); err3=vecnorm([ux3;uy3;uz3]-Uex3)/max(vecnorm(Uex3));
  inb=d3<0.1; errmax(kk)=max(err3);
  fprintf('Stokes SLP all-modes BVP (delta mex): N=%d, np=%d, pmodes=%d, M3=%d\n',N,np,pmodes,M3);
  if any(inb), fprintf('  near band (d<0.1):  max err = %.3e  (%d pts)\n',max(err3(inb)),nnz(inb)); end
  fprintf('  far field (d>=0.1): max err = %.3e  (%d pts)\n',max(err3(~inb)),nnz(~inb));
end

figure(1),clf; semilogy(np_vals,errmax,'o-k'); xlabel('n_p'); ylabel('max err'); grid on;
title('Stokes SLP all-modes BVP (delta mex): h-refinement, p_{modes}=2 n_p');

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
xlabel('x'); ylabel('y'); zlabel('z'); title('Stokes SLP all-modes: 3D target-grid velocity error');

% exportgraphics(figure(1),'axissymsstok_stok_slp_convergence.png','Resolution',200)
% exportgraphics(figure(2),'axissymsstok_stok_slp_error.png','Resolution',200)

function U=stk(X,y,F)
d=X-y; r=vecnorm(d); U=(1/(8*pi))*(F./r+(F.'*d).*d./r.^3);
end
