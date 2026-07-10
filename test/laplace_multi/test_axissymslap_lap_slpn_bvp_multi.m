clearvars; format short e;
addpath('../../utils');
addpath('../../matlab');
addpath('../../external/fmm3d/matlab');

p=16; np_vals=2:10; errmax=nan(1,numel(np_vals));
rot=@(u,th) cos(th)*eye(3)+sin(th)*[0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0]+(1-cos(th))*(u*u');
Zf={@(t) 0.5*sin(t)-1i*1.0*cos(t), @(t) 1.0*sin(t)-1i*0.5*cos(t)}; ac=[0.5 1.0; 1.0 0.5];
yloc={[0.10 0.05; 0 0.05; 0.20 -0.30],[0.20 -0.10; 0 0.15; 0.10 -0.10]}; qq=[1.0 -0.7 0.6 -0.5];

for kk=1:numel(np_vals)
  np=np_vals(kk); M=2*np; nphi=2*M+1; ph=2*pi*(0:nphi-1)/nphi; iside=1; iclosed=0;
  chi_crit=cosh(log(1e13)/M);
  rng(1);

  % 2 particles, random pose; 3D node positions X{k} and 3D node normals NX{k}
  C=cell(2,1); R=cell(2,1); X=cell(2,1); NX=cell(2,1); sq=cell(2,1); geo=cell(2,1);
  for k=1:2
    s=[]; s.p=p; s.Z=Zf{k}; s.tpan=linspace(0,pi,np+1)'; s=quadr(s,[],'p','G');
    d=randn(3,1); C{k}=(k-1)*1.6*d/norm(d); uu=randn(3,1); uu=uu/norm(uu); R{k}=rot(uu,2*pi*rand);
    rho=real(s.x); z=imag(s.x); nr=real(s.nx); nz=imag(s.nx);
    X{k}=R{k}*[(kron(cos(ph(:)),rho)).';(kron(sin(ph(:)),rho)).';(kron(ones(nphi,1),z)).']+C{k};
    NX{k}=R{k}*[(kron(cos(ph(:)),nr)).';(kron(sin(ph(:)),nr)).';(kron(ones(nphi,1),nz)).'];   % 3D node normals
    sq{k}=axisym_to_3d_quadrature(real(s.x),imag(s.x),s.ws,nphi);
    g=[]; g.sx=s.x(:); g.snx=s.nx(:); g.sws=s.ws(:); g.swxp=s.wxp(:); g.tpan=s.tpan(:); g.sxlo=s.xlo(:); g.sxhi=s.xhi(:);
    parc=0; for kj=1:np, parc=max(parc,sum(g.sws((kj-1)*p+(1:p)))); end
    g.rnear=max(max(real(g.sx))*sqrt(2*(chi_crit-1)), 2*parc); geo{k}=g;
  end
  Nmer=numel(geo{1}.sx); N=Nmer*nphi;
  Xall=[X{1} X{2}]; Nall=[NX{1} NX{2}]; wall=[sq{1}.w sq{2}.w]; Mtot=size(Xall,2); blk={1:N, N+1:2*N};

  % interior point charges -> exact exterior potential uex (eval) and exact surface flux dn(uex) (Neumann data)
  Y=zeros(3,0); for k=1:2, Y=[Y, R{k}*yloc{k}+C{k}]; end
  uex=@(P) qq(1)./(4*pi*vecnorm(P-Y(:,1)))+qq(2)./(4*pi*vecnorm(P-Y(:,2))) ...
          +qq(3)./(4*pi*vecnorm(P-Y(:,3)))+qq(4)./(4*pi*vecnorm(P-Y(:,4)));
  gflux=zeros(Mtot,1);                                            % dn(uex) = -sum_k q_k (n.(X-Y_k))/(4pi r^3)
  for k=1:numel(qq), dXY=Xall-Y(:,k); r=vecnorm(dXY); gflux=gflux - qq(k)*(sum(Nall.*dXY,1).')./(4*pi*r(:).^3); end

  % --- (1) naive FULL S' (normal-derivative of SLP) operator (FMM-replaceable; diag 0) ---
  S=lapStrac(Xall,Nall,Xall,wall); S(1:Mtot+1:end)=0;
  Scorr=sparse(Mtot,Mtot);

  % --- (2) self DIFFERENCE correction: (accurate S' self) - (naive S' self, diag 0) ---
  for k=1:2
    g=geo{k};
    [Ab,F,Finv,~,~]=axp_physmat_mex(1,2,1,p,np,g.sx,g.snx,g.sws,g.swxp,g.tpan,g.sxlo,g.sxhi,M,iside,iclosed,1.0,eye(3));
    Sns=lapStrac(X{k},NX{k},X{k},sq{k}.w); Sns(1:N+1:end)=0;
    Scorr(blk{k},blk{k})=Scorr(blk{k},blk{k})+(real(Finv*Ab*F)-Sns);
  end

  % --- (3) off-diagonal near DIFFERENCE correction (kdtree): SLPn flux close - naive ---
  for pq=1:2
    gs=geo{pq}; other=setdiff(1:Mtot, blk{pq});
    nb=~cellfun(@isempty, rangesearch(X{pq}', Xall(:,other)', gs.rnear));
    cand=other(nb); if isempty(cand), continue; end
    Pn=Xall(:,cand); phir=atan2(Pn(2,:),Pn(1,:)); txn=(hypot(Pn(1,:),Pn(2,:))+1i*Pn(3,:)).'; tphin=phir.';
    NPc=Nall(:,cand); nrho=NPc(1,:).*cos(phir)+NPc(2,:).*sin(phir); nphic=-NPc(1,:).*sin(phir)+NPc(2,:).*cos(phir);
    tnxc=(nrho+1i*NPc(3,:)).'; tnphic=nphic.';
    [B,Fs,nmask]=axp_offdiagphysmat_mex(1,2,1,numel(cand),txn,tphin,[0;0;0],eye(3),p,np,gs.sx,gs.snx,gs.sws,gs.swxp,gs.tpan,gs.sxlo,gs.sxhi,C{pq},R{pq},M,iside,iclosed,1.0,0,tnxc,tnphic);
    nr=cand(nmask); Bf=real(B*Fs);
    Sno=lapStrac(Xall(:,nr),Nall(:,nr),X{pq},sq{pq}.w);
    Scorr(nr,blk{pq})=Scorr(nr,blk{pq})+(Bf(nmask,:)-Sno);
  end

  % --- Neumann flux data + solve ---
  sigma=(S+Scorr)\gflux; sig={sigma(blk{1}), sigma(blk{2})};

  % --- 3D exterior target grid (outside both particles) ---
  cmid=(C{1}+C{2})/2; gv=linspace(-3,3,22); [GX,GY,GZ]=ndgrid(gv,linspace(-1,3,12),gv); P=[GX(:)';GY(:)';GZ(:)']+cmid;
  ext=true(1,size(P,2)); dmin=inf(1,size(P,2));
  for k=1:2, loc=R{k}.'*(P-C{k}); sd=(hypot(loc(1,:),loc(2,:))/ac(k,1)).^2+(loc(3,:)/ac(k,2)).^2; ext=ext&(sd>1.05); dmin=min(dmin,sd-1); end
  Pe=P(:,ext); M3=size(Pe,2); dmin=dmin(ext);

  % --- field eval: u = S[sigma] -> naive Lap3dSLPmat + per-source near SLP-value close correction ---
  u=Lap3dSLPmat(struct('x',Pe), struct('x',Xall,'w',wall))*real(sigma);
  for pq=1:2
    gs=geo{pq};
    cand=find(~cellfun(@isempty, rangesearch(X{pq}', Pe', gs.rnear))); if isempty(cand), continue; end
    Pn=Pe(:,cand); txn=(hypot(Pn(1,:),Pn(2,:))+1i*Pn(3,:)).'; tphin=atan2(Pn(2,:),Pn(1,:)).';
    [B,Fs,nmask]=axp_offdiagphysmat_mex(1,1,1,numel(cand),txn,tphin,[0;0;0],eye(3),p,np,gs.sx,gs.snx,gs.sws,gs.swxp,gs.tpan,gs.sxlo,gs.sxhi,C{pq},R{pq},M,iside,iclosed,1.0,0);
    unaive=Lap3dSLPmat(struct('x',Pn), struct('x',X{pq},'w',sq{pq}.w))*real(sig{pq});
    uclose=real(B*(Fs*sig{pq}));
    ci=cand(nmask); u(ci)=u(ci)-unaive(nmask)+uclose(nmask);
  end

  Uex=uex(Pe).'; err=abs(u-Uex)/max(abs(Uex)); inb=dmin.'<0.1; errmax(kk)=max(err);
  fprintf('np=%2d  N=%5d/part  M=%2d  M3=%5d :  near(d<0.1) %.2e   far(d>=0.1) %.2e\n', ...
          np, N, M, M3, max([err(inb);0]), max(err(~inb)));

  figure(2),clf; hold on;
  scatter3(Pe(1,:),Pe(2,:),Pe(3,:),14,log10(max(err,1e-17)),'filled');
  plot3(X{1}(1,:),X{1}(2,:),X{1}(3,:),'.','Color',[.6 .6 .6]);
  plot3(X{2}(1,:),X{2}(2,:),X{2}(3,:),'.','Color',[.6 .6 .6]);
  plot3(Y(1,:),Y(2,:),Y(3,:),'p','MarkerSize',16,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
  cb=colorbar; cb.Label.String='log_{10}|u_h-u_{exact}|/max|u_{exact}|';
  clim([-16 -4]); axis equal; view(35,18); grid on; colormap('jet');
  xlabel('x'); ylabel('y'); zlabel('z'); title('multi-particle Laplace SLPn (S'', close-correction): 3D field error (last refinement)');
end

% --- convergence plot ---
figure(1),clf; semilogy(np_vals,errmax,'o-k'); grid on; xlabel('n_p'); ylabel('max field err');
title('multi-particle Laplace SLPn (S'', close-correction): h-refinement, p=16, M=2 n_p');

function T = lapStrac(Xt, Nt, xs, ws)
% S' (Laplace single-layer normal-derivative) operator: target n-deriv matrix.  Needs target normals Nt.
[~,An] = Lap3dSLPmat(struct('x',Xt,'nx',Nt), struct('x',xs,'w',ws)); T = real(An);
end
