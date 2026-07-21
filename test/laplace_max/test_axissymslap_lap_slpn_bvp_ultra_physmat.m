clearvars; format short e;
addpath('../../../axibie/utils');
addpath('../../utils');
addpath('../../matlab');
addpath('../../external/fmm3d/matlab');                       % lfmm3d / Lap3dSLPfmm

% ULTRA single case (K=15^3=3375, gap 0.4) of the Laplace SLPn (S') exterior-Neumann
% benchmark on the handle framework.  SCALAR (nc=1): no mu, no R sandwich, no deflation
% ((-1/2 + S') is full-rank second-kind); ~1/9 the Stokes footprint.

global APPLY_T APPLY_N

p=16; np=4; gate=2.0; iside=1; iclosed=0;
M=2*np; nphi=2*M+1;
Kside=15; K=Kside^3;
fmmeps=1e-12;
rot=@(u,th) cos(th)*eye(3)+sin(th)*[0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0]+(1-cos(th))*(u*u');
Zf={@(t) 0.5*sin(t)-1i*1.0*cos(t), @(t) 1.0*sin(t)-1i*0.5*cos(t)};
ac=[0.5 1.0; 1.0 0.5];
Rbound=1.0;
yloc={[0.10 0.05; 0 0.05; 0.20 -0.30],[0.20 -0.10; 0 0.15; 0.10 -0.10]};
qloc={[1.0 -0.7],[0.8 -0.6]};

geo=cell(2,1); sq=cell(2,1);
for sh=1:2
  s=[]; s.p=p; s.Z=Zf{sh}; s.tpan=linspace(0,pi,np+1)'; s=quadr(s,[],'p','G');
  g=[]; g.sx=s.x(:); g.snx=s.nx(:); g.sws=s.ws(:); g.swxp=s.wxp(:); g.tpan=s.tpan(:);
  g.rho=real(s.x); g.z=imag(s.x); g.nr=real(s.nx); g.nz=imag(s.nx); geo{sh}=g;
  sq{sh}=axisym_to_3d_quadrature(g.rho,g.z,s.ws,g.nr,g.nz,nphi);
end
N=numel(geo{1}.sx); Nnod=N*nphi;
rnear=gate*sqrt(max(arrayfun(@(j) sum(geo{1}.sws((j-1)*p+(1:p))), 1:np)));
Lsp=2*Rbound + 0.4;
ncores=feature('numcores');
fprintf('ULTRA2: K=%d  p=%d np=%d nphi=%d (Nnod=%d/particle, dof=%d)  ncores=%d\n', ...
        K, p, np, nphi, Nnod, K*Nnod, ncores);
rng(1);

% 1. lattice: alternating shapes, random rotations
C=cell(K,1); R=cell(K,1); ty=zeros(K,1); X=cell(K,1); NX=cell(K,1);
k=0;
for i1=1:Kside, for i2=1:Kside, for i3=1:Kside
  k=k+1; ty(k)=mod(k,2)+1;
  C{k}=Lsp*([i1;i2;i3]-(Kside+1)/2);
  uu=randn(3,1); uu=uu/norm(uu); R{k}=rot(uu,2*pi*rand);
  sqk=sq{ty(k)};
  X{k}=R{k}*sqk.x+C{k};
  NX{k}=R{k}*sqk.nx;
end, end, end
Xall=zeros(3,K*Nnod); Nall=Xall; wall=zeros(1,K*Nnod);
blk=cell(K,1);
for k=1:K
  nodeidx=(k-1)*Nnod+(1:Nnod); blk{k}=nodeidx;
  Xall(:,nodeidx)=X{k}; Nall(:,nodeidx)=NX{k}; wall(nodeidx)=sq{ty(k)}.w;
end

% 2. interior point charges -> exact exterior potential + surface flux
Y=zeros(3,2*K); Q=zeros(1,2*K);
for k=1:K, Y(:,2*k-1:2*k)=R{k}*yloc{ty(k)}+C{k}; Q(2*k-1:2*k)=qloc{ty(k)}; end
uexf=@(P) lapsum(P,Y,Q);
g_tr=zeros(1,K*Nnod);
for j=1:2*K, d=Xall-Y(:,j); r=vecnorm(d); g_tr=g_tr - Q(j)*(sum(Nall.*d,1))./(4*pi*r.^3); end
g_tr=g_tr(:);

% flat per-particle geometry (consumed by every corr_setup below)
sxs=zeros(N,K); snxs=sxs; swss=zeros(N,K); swxps=zeros(N,K); tpans=zeros(np+1,K);
for k=1:K
  gs=geo{ty(k)};
  sxs(:,k)=gs.sx; snxs(:,k)=gs.snx; swss(:,k)=gs.sws; swxps(:,k)=gs.swxp; tpans(:,k)=gs.tpan;
end
pv=p*ones(K,1); npv=np*ones(K,1); pmv=M*ones(K,1);
geomoff=1+(0:K)'*N; tpanoff=1+(0:K)'*(np+1); targoff=1+(0:K)'*Nnod;
nsx=K*N; ntpan=K*(np+1);
sxf=sxs(:); snxf=snxs(:); swsf=swss(:); swxpf=swxps(:); tpanf=tpans(:);

% 3. self corrections (laplace SLPn, iinter=1)
tself=tic;
[hself, nbself] = axpso_corr_setup_mex(1, 2, 1, 0, 0.0, 1, K, pv, npv, pmv, iside, iclosed, gate, ...
          geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
          Rbound+rnear, K*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);
tself=toc(tself); selfMB=nbself/1e6;

% 4. cross corrections (iinter=2)
tcross=tic;
[hcross, nbcross] = axpso_corr_setup_mex(1, 2, 1, 0, 0.0, 2, K, pv, npv, pmv, iside, iclosed, gate, ...
           geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
           Rbound+rnear, K*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);
tcross=toc(tcross); crossMB=nbcross/1e6;

% 5. matrix-free operator (no deflation: second-kind full-rank)
APPLY_T=0; APPLY_N=0;
applyA=@(x) applytimed(x, fmmeps, Xall, Nall, wall, hself, hcross);

% 6. solve
tso=tic;
[sigma,flag,relres,iter]=gmres(applyA, g_tr, [], 1e-9, 200);
tsolve=toc(tso); Teval=APPLY_T/APPLY_N;

% 7. off-surface eval: near/mid/far check-point shells + field-eval handle (ilayer=1)
teo=tic;
rng(2); nchk=min(200,Nnod); Pe=[];
for d=[0.05 0.3 1.0], isel=randperm(K*Nnod,nchk); Pe=[Pe, Xall(:,isel)+d*Nall(:,isel)]; end %#ok<AGROW>
keep=true(1,size(Pe,2));
for k=1:K
  loc=R{k}.'*(Pe-C{k});
  sd=(hypot(loc(1,:),loc(2,:))/ac(ty(k),1)).^2+(loc(3,:)/ac(ty(k),2)).^2;
  keep=keep & (sd>1.10);
end
Pe=Pe(:,keep); Me=size(Pe,2);
heval = axpso_corr_setup_mex(1, 1, 1, 0, 0.0, 2, K, pv, npv, pmv, iside, iclosed, gate, ...
          geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
          Rbound+rnear, Me, ones(K+1,1), Pe, Pe, cat(3,R{:}), [C{:}], 0);
u=Lap3dSLPfmm(struct('x',Pe), struct('x',Xall,'w',wall), sigma, fmmeps);
u=axpso_corr_apply_mex(heval, K*Nnod, sigma, Me, u);
tevaloff=toc(teo);
Ue=uexf(Pe).'; err=abs(u-Ue)/max(abs(Ue));
dsurf=zeros(1,Me); for i=1:Me, dsurf(i)=min(vecnorm(Xall-Pe(:,i))); end
fprintf('  err by shell:  d~0.05: %.2e   d~0.3: %.2e   d~1.0: %.2e\n', ...
        max(err(dsurf<0.15)), max(err(dsurf>=0.15 & dsurf<0.6)), max(err(dsurf>=0.6)));

% ---- visualization: particle SURFACES (one combined coarse patch; the full 3.7M-node
% cloud would choke MATLAB) + field-error scatter at the check points.  Stokeslet
% markers dropped (particles are tiny at this scale).
dsurf=zeros(1,Me);
for i=1:Me, dsurf(i)=min(vecnorm(Xall-Pe(:,i))); end          % dist to nearest surface node
msel=1:4:N; nm=numel(msel);                                   % meridian decimation for display
F0=zeros((nm-1)*nphi,4); c=0;                                 % one particle's quad faces
for a=1:nphi
  a2=mod(a,nphi)+1;                                           % azimuthal wrap
  for j=1:nm-1
    c=c+1; F0(c,:)=[(a-1)*nm+j, (a-1)*nm+j+1, (a2-1)*nm+j+1, (a2-1)*nm+j];
  end
end
VV=zeros(3,K*nm*nphi);
for k=1:K
  VV(:,(k-1)*nm*nphi+(1:nm*nphi))=X{k}(:,reshape(msel(:)+N*(0:nphi-1),1,[]));
end
FF=repmat(F0,K,1)+repelem((0:K-1)'*(nm*nphi), size(F0,1))*ones(1,4);
figure(2),clf; hold on;
patch('Vertices',VV.','Faces',FF,'FaceColor',[0.75 0.75 0.78],'EdgeColor','none');
scatter3(Pe(1,:),Pe(2,:),Pe(3,:),14,log10(max(err,1e-17)),'filled');
cb=colorbar; cb.Label.String='log_{10}|u_h-u_{exact}|/max|u_{exact}|';
axis equal; view(35,18); grid on; colormap('jet');
camlight headlight; lighting gouraud; material dull;
title(sprintf('K=%d particle surfaces + field error at check points',K)); drawnow;
% exportgraphics(figure(2),sprintf('axissymslap_lap_slpn_ultra_physmat_K%d.png',K),'Resolution',200)
fprintf('  err by shell:  d~0.05: %.2e   d~0.3: %.2e   d~1.0: %.2e\n', ...
        max(err(dsurf<0.15)), max(err(dsurf>=0.15 & dsurf<0.6)), max(err(dsurf>=0.6)));

% 8. the data point
dof=K*Nnod;
fprintf(['K=%4d  dof=%8d :  T_setup %6.1fs (T_self %6.1fs %.0f MB + T_cross %6.1fs %.0f MB)', ...
         '   [%d threads]   %8.0f dof/s/core\n'], ...
        K, dof, tself+tcross, tself, selfMB, tcross, crossMB, ncores, dof/(tself+tcross)/ncores);
fprintf(['                        T_eval %6.3f s/iter (fmm eps %.0e, %4.0f us/src)  iters %3d (flag=%d, relres=%.1e)  ', ...
         'T_solve %5.0fs  T_eval_off %5.1fs (%d tgts)  err %.2e   [%d cores]  %8.0f dof/s/core\n'], ...
        Teval, fmmeps, Teval/(K*Nnod)*1e6, iter(2), flag, relres, tsolve, tevaloff, Me, max(err), ...
        ncores, dof/Teval/ncores);

function y=applytimed(x, fmmeps, Xall, Nall, wall, hself, hcross)
global APPLY_T APPLY_N
t0=tic;
srcinfo.sources=Xall; srcinfo.charges=(x(:).').*wall;
U=lfmm3d(fmmeps,srcinfo,2);
y=(sum(U.grad.*Nall,1)).';
n=numel(x);
y=axpso_corr_apply_mex(hself,  n, x, n, y);
y=axpso_corr_apply_mex(hcross, n, x, n, y);
APPLY_T=APPLY_T+toc(t0); APPLY_N=APPLY_N+1;
end
function u=lapsum(P,Y,Q)
u=zeros(1,size(P,2));
for j=1:size(Y,2), u=u+Q(j)./(4*pi*vecnorm(P-Y(:,j))); end
end
