clearvars; format short e;
addpath('../../utils');
addpath('../../matlab');
addpath('../../external/fmm3d/matlab');                       % stfmm3d / Sto3dSLPfmm

global APPLY_T APPLY_N                                        % T_eval accumulator (set in applytimed)

p=16; np=4; mu=1.0; gate=2.0; iside=1; iclosed=0;
M=2*np; nphi=2*M+1;
Ksides=[2 3 4];                                                 % K = 8, 27 (K=64/125: user run, per task-5 brief)
rot=@(u,th) cos(th)*eye(3)+sin(th)*[0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0]+(1-cos(th))*(u*u');
Zf={@(t) 0.5*sin(t)-1i*1.0*cos(t), @(t) 1.0*sin(t)-1i*0.5*cos(t)};   % 2 particle shapes (prolate/oblate)
ac=[0.5 1.0; 1.0 0.5];                                        % semi-axes (rho, z) per shape
Rbound=1.0;                                                   % bounding radius of both shapes
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
Lsp=2*Rbound + 0.6;                                           % lattice: worst-case surface gap 0.6
fmmeps=1e-12;                                                 % FMM tolerance (matvec + off-surface eval)
ncores=feature('numcores');                                   % physical cores (dof/s/core normalization)
fprintf('p=%d np=%d nphi=%d (Nnod=%d/particle)  gate=%.1f  rnear=%.2f  lattice spacing %.2f  ncores=%d\n', ...
        p, np, nphi, Nnod, gate, rnear, Lsp, ncores);

for Kside=Ksides
  K=Kside^3;
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
  g_tr=g_tr(:);                                               % exact traction, node-interleaved (lab)

  % 3. SELF corrections (NEW handle framework): stacked (N,K) meridian geometry -> flat
  %    per-particle layout (geomoff/tpanoff), ONE handle-based corr_setup call.  ikernel=2
  %    (stokes), ilayer=2 (SLPn), params=mu (real part read as viscosity), iinter=1 (self).
  %    closesize + closestokslpn per particle happen INSIDE the module, addressed by hself.
  sxs=zeros(N,K); snxs=sxs; swss=zeros(N,K); swxps=zeros(N,K); tpans=zeros(np+1,K);
  for k=1:K
    gs=geo{ty(k)};
    sxs(:,k)=gs.sx; snxs(:,k)=gs.snx; swss(:,k)=gs.sws; swxps(:,k)=gs.swxp; tpans(:,k)=gs.tpan;
  end
  pv=p*ones(K,1); npv=np*ones(K,1); pmv=M*ones(K,1);
  geomoff=1+(0:K)'*N; tpanoff=1+(0:K)'*(np+1); targoff=1+(0:K)'*Nnod;
  nsx=K*N; ntpan=K*(np+1);
  sxf=sxs(:); snxf=snxs(:); swsf=swss(:); swxpf=swxps(:); tpanf=tpans(:);
  tself=tic;
  [hself, nbself] = axpso_corr_setup_mex(2, 2, mu, 1, K, pv, npv, pmv, iside, iclosed, gate, ...
            geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
            Rbound+rnear, K*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);
  tself=toc(tself); selfMB=nbself/1e6;

  % 4. CROSS corrections (NEW handle framework): ONE handle-based corr_setup call (iinter=2);
  %    the module builds its own canonical kdtree, ball query / own-block drop / rotation into
  %    each source frame / sizing / canonical compaction, addressed by hcross.  Reuses the
  %    flat per-particle geometry built in section 3.
  tcross=tic;
  [hcross, nbcross] = axpso_corr_setup_mex(2, 2, mu, 2, K, pv, npv, pmv, iside, iclosed, gate, ...
             geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
             Rbound+rnear, K*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);
  tcross=toc(tcross); crossMB=nbcross/1e6;

  % 5. matrix-free operator (timed per apply via the wrapper) + rank-K deflation.  The near
  %    correction now rides the handles hself/hcross (the R sandwich is INSIDE the module).
  Nh=zeros(3*K*Nnod,K);
  for k=1:K, Nh(blk{k},k)=NX{k}(:); Nh(:,k)=Nh(:,k)/norm(Nh(:,k)); end
  APPLY_T=0; APPLY_N=0;
  applyA=@(x) applytimed(x, Nh, fmmeps, Xall, Nall, wall, hself, hcross);

  % 6. solve;  T_eval = mean apply time over ALL applies inside the solve
  tso=tic;
  [sigma,flag,relres,iter]=gmres(applyA, g_tr, [], 1e-9, 200);
  tsolve=toc(tso);
  Teval=APPLY_T/APPLY_N;

  % 7. off-surface TARGET SETUP: regular ng^3 grid, exact spheroid inside-filter
  %    (replaces physmat2's random check-point shells).  The evaluation itself --
  %    Sto3dSLPfmm seed + ONE fortran fieldeval call -- is added next.
  teo=tic;
  ng=16; margin=0.9;
  glo=-Lsp*(Kside-1)/2 - Rbound - margin;
  xs=linspace(glo,-glo,ng);
  [xg,yg,zg]=ndgrid(xs,xs,xs);                                % first index fastest = x
  Pg=[xg(:).'; yg(:).'; zg(:).'];
  inside=false(1,size(Pg,2));                                 % exact spheroid inside-test
  for k=1:K
    loc=R{k}.'*(Pg-C{k});
    sd=(hypot(loc(1,:),loc(2,:))/ac(ty(k),1)).^2+(loc(3,:)/ac(ty(k),2)).^2;
    inside=inside | (sd<1.0);
  end
  Pe=Pg(:,~inside); Me=size(Pe,2);
  fprintf('  targets: %d^3 grid, %d outside particles (h=%.3f)\n', ng, Me, xs(2)-xs(1));
  % exact solution
  Ue = zeros(3,size(Pg,2));
  Ue(:,~inside) = uexf(Pg(:,~inside));
  sigblk=reshape(reshape(sigma,3,[]).',[],1);                 % node-interleaved -> component-blocked
  ublk=Sto3dSLPfmm(struct('x',Pe), struct('x',Xall,'w',wall), sigblk, fmmeps);  % FMM velocity at outside targets
  Uh=zeros(3,size(Pg,2)); Uh(:,~inside)=reshape(ublk,[],3).'; % full-grid shape (inside zeroed)
  uc=reshape(Uh(:,~inside),[],1);                             % seed with naive velocity (node-interleaved)
  % FIELD-EVAL correction (NEW handle framework): reuse the cross path (iinter=2) with ilayer=1
  % (SLP VELOCITY) and targets=Pe; targoff=ones(K+1,1) -> NO own-block drop.  ONE apply, source
  % 3*K*Nnod -> target 3*Me (the lab<->local R sandwich is inside the module).
  heval = axpso_corr_setup_mex(2, 1, mu, 2, K, pv, npv, pmv, iside, iclosed, gate, ...
            geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
            Rbound+rnear, Me, ones(K+1,1), Pe, Pe, cat(3,R{:}), [C{:}], 0);
  uc=axpso_corr_apply_mex(heval, 3*K*Nnod, sigma, 3*Me, uc);
  Uh(:,~inside)=reshape(uc,3,[]);
  tevaloff=toc(teo);
  Ue=uexf(Pe);
  scal=max(vecnorm(Ue));
  errg=zeros(1,size(Pg,2)); errg(~inside)=vecnorm(Uh(:,~inside)-Ue)/scal;   % per-point (step 9 field)
  err=max(errg);

  % 8. the data point (paper-style, two lines: SETUP is the OpenMP solver-module pass
  %    (mex built with OMP=ON -- default libomp thread count = ncores); EVAL/SOLVE is
  %    dominated by the multithreaded FMM.  dof = 3 x quadrature points (Stokes
  %    3-vector); dof/s/core normalizes each phase by the cores it actually uses.
  dof=3*K*Nnod;
  fprintf(['K=%4d  dof=%8d :  T_setup %6.1fs (T_self %6.1fs %.0f MB + T_cross %6.1fs %.0f MB)', ...
           '   [%d threads]   %8.0f dof/s/core\n'], ...
          K, dof, tself+tcross, tself, selfMB, tcross, crossMB, ncores, dof/(tself+tcross)/ncores);
  fprintf(['                        T_eval %6.3f s/iter (fmm eps %.0e, %4.0f us/src)  iters %3d (flag=%d, relres=%.1e)  ', ...
           'T_solve %5.0fs  T_eval_off %5.1fs (%d tgts)  err %.2e   [%d cores]  %8.0f dof/s/core\n'], ...
          Teval, fmmeps, Teval/(K*Nnod)*1e6, iter(2), flag, relres, tsolve, tevaloff, Me, max(err), ...
          ncores, dof/Teval/ncores);

end

function y=applytimed(x, Nh, fmmeps, Xall, Nall, wall, hself, hcross)
% timed K-particle matrix-free (-1/2 I + S') matvec: global stfmm3d naive traction (self term
% excluded) + the handle-based near correction (self then cross, each ACCUMULATING; the lab<->
% local R sandwich happens INSIDE the module) + rank-K pressure-gauge deflation Nh*(Nh'*x).
% T_eval = accumulated apply time / count (per-iteration average, excludes gmres orthog).
global APPLY_T APPLY_N
t0=tic;
srcinfo.sources=Xall; srcinfo.stoklet=reshape(x,3,[]).*wall;
U=stfmm3d(fmmeps,srcinfo,3);
Gs=U.grad+permute(U.grad,[2 1 3]);
y=squeeze(sum(Gs.*reshape(Nall,1,3,[]),2)) - U.pre.*Nall;
y=y(:);
n=numel(x);
y=axpso_corr_apply_mex(hself,  n, x, n, y);                  % self  correction (accumulate; nu==nx)
y=axpso_corr_apply_mex(hcross, n, x, n, y);                  % cross correction (accumulate; nu==nx)
y=y + Nh*(Nh'*x);                                            % rank-K pressure-gauge deflation
APPLY_T=APPLY_T+toc(t0); APPLY_N=APPLY_N+1;
end
function u=stksum(P,Y,Fg,stk)
u=zeros(3,size(P,2));
for j=1:size(Y,2), u=u+stk(P,Y(:,j),Fg(:,j)); end
end
function T=stktrac(X,n,y,F)
d=X-y; r=vecnorm(d); rF=F.'*d; rn=sum(n.*d,1); T=-(3/(4*pi))*(rF.*rn).*d./r.^5;
end
