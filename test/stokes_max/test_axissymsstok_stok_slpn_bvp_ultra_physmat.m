clearvars; format short e;
addpath('../../../axibie/utils');
addpath('../../utils');
addpath('../../matlab');
addpath('../../external/fmm3d/matlab');                       % stfmm3d / Sto3dSLPfmm

% ULTRA single-case run of the stokes_max SLPn benchmark: K = 15^3 = 3375 particles,
% np=4 (Nnod=1088/particle, ~3.7 M surface nodes, ~11 M dof).  Sized for a 512 GB
% machine (Mac Studio M3 Ultra): MEASURED 2026-07-03 at gap 0.4 (RESULTS_ULTRA.md)
% cross 73.4 + self 4.56 MB/particle -> ~263 GB corrections; DO NOT run on 64 GB.
% SOLVER-MODULE pipeline (test_axissymsstok_stok_slpn_bvp_max_physmat2.m conventions):
% the three K-loops are single axpso_ mex calls with the correction blocks resident in
% Fortran module state (never crossing the MATLAB boundary); build the mex with
% `make mex OMP=ON` for the OpenMP setup/apply (serial otherwise).  Scale adaptation
% kept: SPARSE rank-K deflation matrix (dense would be ~300 GB).
% LIFETIME RULE: do NOT `clear mex` / `clear all` between the setup calls and the
% solve/eval -- that unloads the gateway and destroys the module-resident blocks
% (axpso_stokslpn_release_mex(0) at the end frees them explicitly).

global APPLY_T APPLY_N                                        % T_eval accumulator (set in applytimed)

p=16; np=4; mu=1.0; gate=2.0; iside=1; iclosed=0;
M=2*np; nphi=2*M+1;
Kside=15; K=Kside^3;                                          % the single ultra case
fmmeps=1e-12;                                                 % FMM tolerance (matvec + off-surface eval)
rot=@(u,th) cos(th)*eye(3)+sin(th)*[0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0]+(1-cos(th))*(u*u');
Zf={@(t) 0.5*sin(t)-1i*1.0*cos(t), @(t) 1.0*sin(t)-1i*0.5*cos(t)};   % 2 particle shapes (prolate/oblate)
ac=[0.5 1.0; 1.0 0.5];                                        % semi-axes (rho, z) per shape
Rbound=1.0;
stk=@(P,y,Fv)(1/(8*pi))*( Fv./vecnorm(P-y) + (Fv.'*(P-y)).*(P-y)./vecnorm(P-y).^3 );
yloc={[0.10 0.05; 0 0.05; 0.20 -0.30],[0.20 -0.10; 0 0.15; 0.10 -0.10]};
Floc={[1.0 -0.5; -0.7 0.6; 0.5 -0.4],[0.6 -0.4; 0.3 -0.5; -0.7 0.5]};

% ---- per-shape meridian geometry (construction only -- self passes are PER PARTICLE) ----
geo=cell(2,1); sq=cell(2,1);
for sh=1:2
  s=[]; s.p=p; s.Z=Zf{sh}; s.tpan=linspace(0,pi,np+1)'; s=quadr(s,[],'p','G');
  g=[]; g.sx=s.x(:); g.snx=s.nx(:); g.sws=s.ws(:); g.swxp=s.wxp(:); g.tpan=s.tpan(:);
  g.rho=real(s.x); g.z=imag(s.x); g.nr=real(s.nx); g.nz=imag(s.nx); geo{sh}=g;
  sq{sh}=axisym_to_3d_quadrature(g.rho,g.z,s.ws,g.nr,g.nz,nphi);
end
N=numel(geo{1}.sx); Nnod=N*nphi;
rnear=gate*sqrt(max(arrayfun(@(j) sum(geo{1}.sws((j-1)*p+(1:p))), 1:np)));  % near-zone width
Lsp=2*Rbound + 0.4;                                           % lattice: worst-case surface gap 0.4 (tighter than
                                                              % the K<=125 runs' 0.6; still >> lubrication regime)
ncores=feature('numcores');
fprintf('ULTRA: K=%d  p=%d np=%d nphi=%d (Nnod=%d/particle, dof=%d)  ncores=%d\n', ...
        K, p, np, nphi, Nnod, 3*K*Nnod, ncores);
fprintf('predicted (measured gap-0.4 constants): corrections ~%.0f GB, setup ~%.0f min serial / ~%.0f min at %d threads\n', ...
        K*(73.4+4.56)/1024, K*0.577/60, K*0.577/60/ncores, ncores);
rng(1);

% 1. K particles on the lattice: alternating shapes, random rotations
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
  nodeidx=(k-1)*Nnod+(1:Nnod); blk{k}=(k-1)*3*Nnod+(1:3*Nnod);
  Xall(:,nodeidx)=X{k}; Nall(:,nodeidx)=NX{k}; wall(nodeidx)=sq{ty(k)}.w;
end

% 2. interior Stokeslets (2 per particle) -> exact velocity + surface traction data
Y=zeros(3,2*K); Fg=zeros(3,2*K);
for k=1:K, Y(:,2*k-1:2*k)=R{k}*yloc{ty(k)}+C{k}; Fg(:,2*k-1:2*k)=R{k}*Floc{ty(k)}; end
uexf=@(P) stksum(P,Y,Fg,stk);
g_tr=zeros(3,K*Nnod);
for j=1:2*K, g_tr=g_tr+stktrac(Xall,Nall,Y(:,j),Fg(:,j)); end
g_tr=g_tr(:);                                                 % exact traction, node-interleaved (lab)

% 3. SELF corrections: stacked (N,K) geometry, ONE module-resident selfsetup mex call
%    (blocks stored as Fortran module state; OpenMP over particles with OMP=ON).
tself=tic;
sxs=zeros(N,K); snxs=sxs; swss=zeros(N,K); swxps=zeros(N,K); tpans=zeros(np+1,K);
for k=1:K
  gs=geo{ty(k)};
  sxs(:,k)=gs.sx; snxs(:,k)=gs.snx; swss(:,k)=gs.sws; swxps(:,k)=gs.swxp; tpans(:,k)=gs.tpan;
end
[ntcxs,nbself]=axpso_stokslpn_selfsetup_mex(K, p, np, nphi, M, iside, iclosed, gate, mu, ...
    sxs, snxs, swss, swxps, tpans, zeros(K,1), 0);
tself=toc(tself); selfMB=nbself/1e6;
fprintf('  selfsetup done: %.0fs, %.1f GB\n', tself, nbself/1e9);

% 4. CROSS corrections: ONE module-resident crosssetup mex call (the module builds its
%    own canonical tree, does ball query / own-block drop / rotation / sizing /
%    canonical mapping / unique compaction, stores compact per-source blocks).
tcross=tic;
[ntcxx,nuo,npair,nbcross]=axpso_stokslpn_crosssetup_mex(K, p, np, nphi, M, iside, iclosed, gate, mu, ...
    Rbound+rnear, Xall, Nall, cat(3,R{:}), [C{:}], sxs, snxs, swss, swxps, tpans, ...
    zeros(K,1), zeros(K,1), 0, 0);
tcross=toc(tcross); crossMB=nbcross/1e6;
Rcat=cat(3,R{:});
fprintf('  crosssetup done: %.0fs, %.1f GB, %d src particles with near pairs\n', tcross, nbcross/1e9, npair);

% 5. matrix-free operator + SPARSE rank-K deflation (dense Nh would be ~300 GB here)
nhv=zeros(3*K*Nnod,1);
for k=1:K, nhv(blk{k})=NX{k}(:)/norm(NX{k}(:)); end
Nh=sparse(1:3*K*Nnod, repelem(1:K,3*Nnod), nhv, 3*K*Nnod, K);
APPLY_T=0; APPLY_N=0;
applyA=@(x) applytimed(x, Nh, fmmeps, Xall, Nall, wall, K, Rcat);

% 6. solve;  T_eval = mean apply time over ALL applies inside the solve
tso=tic;
[sigma,flag,relres,iter]=gmres(applyA, g_tr, [], 1e-9, 60);
tsolve=toc(tso);
Teval=APPLY_T/APPLY_N;

% 7. off-surface evaluation of the FINAL solution (timed)
teo=tic;
rng(2); nchk=200; Pe=[];
for d=[0.05 0.3 1.0]
  isel=randperm(K*Nnod,nchk);
  Pe=[Pe, Xall(:,isel)+d*Nall(:,isel)]; %#ok<AGROW>
end
keep=true(1,size(Pe,2));                                      % exact spheroid inside-test per particle
for k=1:K
  loc=R{k}.'*(Pe-C{k});
  sd=(hypot(loc(1,:),loc(2,:))/ac(ty(k),1)).^2+(loc(3,:)/ac(ty(k),2)).^2;
  keep=keep & (sd>1.10);
end
Pe=Pe(:,keep); Me=size(Pe,2);
sigblk=reshape(reshape(sigma,3,[]).',[],1);
ublk=Sto3dSLPfmm(struct('x',Pe), struct('x',Xall,'w',wall), sigblk, fmmeps);
ucm=reshape(reshape(ublk,[],3).',[],1);
for k=1:K
  gs=geo{ty(k)};
  loc=R{k}.'*(Pe-C{k});
  if min(vecnorm(loc)) > Rbound+rnear, continue, end
  ttx=(hypot(loc(1,:),loc(2,:)) + 1i*loc(3,:)).';
  tcxi=zeros(np+1,1); ntcx=0;
  [tcxi,ntcx]=axps_closesize_mex(Me, ttx, loc, p, np, gs.sx, gs.sws, gate, tcxi, ntcx);
  if ntcx==0, continue, end
  Sk=zeros(3*ntcx,3*nphi*p); ik=zeros(ntcx,1);
  [Sk,ik]=axps_closestokslp_mex(Me, ttx, loc, p, np, nphi, gs.sx, gs.snx, gs.sws, gs.swxp, gs.tpan, gate, ...
      sq{ty(k)}.x, sq{ty(k)}.nx, sq{ty(k)}.w, M, iside, iclosed, mu, ntcx, tcxi, Sk, ik);
  sigk=reshape(R{k}.'*reshape(sigma(blk{k}),3,[]),[],1);
  uc=zeros(3*Me,1);
  uc=axps_closecorrapply_mex(3, Me, p, np, nphi, ntcx, tcxi, ik, Sk, sigk, 0, uc);
  ucm=ucm+reshape(R{k}*reshape(uc,3,[]),[],1);
end
tevaloff=toc(teo);
U=reshape(ucm,3,Me); Ue=uexf(Pe);
err=vecnorm(U-Ue)/max(vecnorm(Ue));

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
% exportgraphics(figure(2),sprintf('axissymsstok_stok_slpn_ultra_physmat_K%d.png',K),'Resolution',200)
fprintf('  err by shell:  d~0.05: %.2e   d~0.3: %.2e   d~1.0: %.2e\n', ...
        max(err(dsurf<0.15)), max(err(dsurf>=0.15 & dsurf<0.6)), max(err(dsurf>=0.6)));

% 8. the data point (paper-style, two lines; T_setup = T_self + T_cross)
dof=3*K*Nnod;
fprintf(['K=%4d  dof=%8d :  T_setup %6.1fs (T_self %6.1fs %.0f MB + T_cross %6.1fs %4d srcs %.0f MB)', ...
         '   [%d threads]   %8.0f dof/s/core\n'], ...
        K, dof, tself+tcross, tself, selfMB, tcross, npair, crossMB, ncores, dof/(tself+tcross)/ncores);
fprintf(['                        T_eval %6.3f s/iter (fmm eps %.0e, %4.0f us/src)  iters %3d (flag=%d, relres=%.1e)  ', ...
         'T_solve %5.0fs  T_eval_off %5.1fs (%d tgts)  err %.2e   [%d cores]  %8.0f dof/s/core\n'], ...
        Teval, fmmeps, Teval/(K*Nnod)*1e6, iter(2), flag, relres, tsolve, tevaloff, Me, max(err), ...
        ncores, dof/Teval/ncores);
ok=axpso_stokslpn_release_mex(0); %#ok<NASGU>                 % free the ~260 GB module state

function y=applytimed(x, Nh, fmmeps, varargin)
% timed operator apply: accumulates total apply time / count for T_eval
global APPLY_T APPLY_N
t0=tic;
y=mvslpnK(x, fmmeps, varargin{:}) + Nh*(Nh'*x);
dt=toc(t0); APPLY_T=APPLY_T+dt; APPLY_N=APPLY_N+1;
fprintf('  gmres apply %3d done: %.1fs (cum %.0fs)\n', APPLY_N, dt, APPLY_T);
end
function u=stksum(P,Y,Fg,stk)
u=zeros(3,size(P,2));
for j=1:size(Y,2), u=u+stk(P,Y(:,j),Fg(:,j)); end
end
function T=stktrac(X,n,y,F)
d=X-y; r=vecnorm(d); rF=F.'*d; rn=sum(n.*d,1); T=-(3/(4*pi))*(rF.*rn).*d./r.^5;
end
function t=mvslpnK(x, fmmeps, Xall, Nall, wall, K, Rcat)
% K-particle matrix-free (-1/2 I + S') matvec: global stfmm3d naive traction (self term
% excluded) + ONE module-resident correction apply (the sigloc build / R sandwich /
% circulant self + compact cross scatter all happen inside the solver module).
srcinfo.sources=Xall; srcinfo.stoklet=reshape(x,3,[]).*wall;
U=stfmm3d(fmmeps,srcinfo,3);
Gs=U.grad+permute(U.grad,[2 1 3]);
t=squeeze(sum(Gs.*reshape(Nall,1,3,[]),2)) - U.pre.*Nall;
t=t(:);
t=axpso_stokslpn_corrapply_mex(K, Rcat, x, t);
end
