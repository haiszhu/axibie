clearvars; format short e;
addpath('../../utils');
addpath('../../matlab');
addpath('../../external/fmm3d/matlab');                       % lfmm3d / Lap3dSLPfmm

global APPLY_T APPLY_N                                        % T_eval accumulator (set in applytimed)

p=16; np=4; gate=2.0; iside=1; iclosed=0;
M=2*np; nphi=2*M+1;
Ksides=[2 3 4 5];                                                 % K = 8, 27, 64, 125
rot=@(u,th) cos(th)*eye(3)+sin(th)*[0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0]+(1-cos(th))*(u*u');
Zf={@(t) 0.5*sin(t)-1i*1.0*cos(t), @(t) 1.0*sin(t)-1i*0.5*cos(t)};   % 2 particle shapes (prolate/oblate)
ac=[0.5 1.0; 1.0 0.5];                                        % semi-axes (rho, z) per shape
Rbound=1.0;                                                   % bounding radius of both shapes
yloc={[0.10 0.05; 0 0.05; 0.20 -0.30],[0.20 -0.10; 0 0.15; 0.10 -0.10]};   % 2 interior charges/particle
qloc={[1.0 -0.7],[0.8 -0.6]};                                % scalar charge strengths per shape

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
fmmeps=1e-12;
ncores=feature('numcores');
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
    nodeidx=(k-1)*Nnod+(1:Nnod); blk{k}=nodeidx;             % SCALAR: one dof per node
    Xall(:,nodeidx)=X{k}; Nall(:,nodeidx)=NX{k}; wall(nodeidx)=sq{ty(k)}.w;
  end

  % for new solver interface
  pv   = p*ones(K,1);  npv = np*ones(K,1);  pmv = M*ones(K,1);
  geomoff = 1 + (0:K)'*N;          % N = p*np   -> flat sx block base per particle
  tpanoff = 1 + (0:K)'*(np+1);     % flat tpan/tcxi block base
  targoff = 1 + (0:K)'*Nnod;       % own-block boundaries (cross own-drop)
  nsx = K*N;  ntpan = K*(np+1);

  % 2. interior point charges (2 per particle) -> exact exterior potential + surface flux data
  Y=zeros(3,2*K); Q=zeros(1,2*K);
  for k=1:K, Y(:,2*k-1:2*k)=R{k}*yloc{ty(k)}+C{k}; Q(2*k-1:2*k)=qloc{ty(k)}; end
  uexf=@(P) lapsum(P,Y,Q);
  g_tr=zeros(1,K*Nnod);                                       % dn(uex) = -sum_j q_j (n.(X-Y_j))/(4pi r^3)
  for j=1:2*K, d=Xall-Y(:,j); r=vecnorm(d); g_tr=g_tr - Q(j)*(sum(Nall.*d,1))./(4*pi*r.^3); end
  g_tr=g_tr(:);                                               % exact surface flux, scalar per node

  % stacked (N,K) meridian geometry -> flat per-particle layout (geomoff/tpanoff), the input
  % the handle-based setup consumes for BOTH self and cross.
  sxs=zeros(N,K); snxs=sxs; swss=zeros(N,K); swxps=zeros(N,K); tpans=zeros(np+1,K);
  for k=1:K
    gs=geo{ty(k)};
    sxs(:,k)=gs.sx; snxs(:,k)=gs.snx; swss(:,k)=gs.sws; swxps(:,k)=gs.swxp; tpans(:,k)=gs.tpan;
  end
  sxf = sxs(:); snxf = snxs(:); swsf = swss(:); swxpf = swxps(:); tpanf = tpans(:);   % (N,K)->flat, (np+1,K)->flat

  % 3. SELF corrections (iinter=1): ONE handle-based setup call -- closesize + closeslpn per
  %    particle happen INSIDE the module, blocks stored as module state, addressed by hself.
  tself=tic;
  hself  = axpso_corr_setup_mex(1, 2, 0.0, 1, K, pv, npv, pmv, iside, iclosed, gate, ...
             geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
             Rbound+rnear, K*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);
  tself=toc(tself);

  % 4. CROSS corrections (iinter=2): ONE handle-based setup call -- the module builds its own
  %    canonical kdtree, ball query / own-block drop / sizing / canonical compaction, stores
  %    the compact per-source blocks as module state, addressed by hcross.
  tcross=tic;
  hcross = axpso_corr_setup_mex(1, 2, 0.0, 2, K, pv, npv, pmv, iside, iclosed, gate, ...
             geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
             Rbound+rnear, K*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);
  tcross=toc(tcross);

  % 5. matrix-free operator (timed per apply).  NO deflation: (-1/2 + S') is second-kind
  %    FULL-RANK for the Laplace exterior-Neumann problem (unlike the Stokes traction
  %    operator, whose per-particle normal null vectors need the rank-K deflator).
  APPLY_T=0; APPLY_N=0;
  applyA=@(x) applytimed(x, fmmeps, Xall, Nall, wall, hself, hcross);

  % 6. solve;  T_eval = mean apply time over ALL applies inside the solve
  tso=tic;
  [sigma,flag,relres,iter]=gmres(applyA, g_tr, [], 1e-9, 200);
  tsolve=toc(tso); Teval=APPLY_T/APPLY_N;

  % 7. off-surface evaluation of the FINAL solution (timed): near/mid/far check points
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
  % FIELD-EVAL correction handle: reuse the cross path (iinter=2) with ilayer=1 (SLP VALUE)
  % and targets=Pe; targoff=ones(K+1,1) makes every own-range empty -> NO own-block drop.
  heval = axpso_corr_setup_mex(1, 1, 0.0, 2, K, pv, npv, pmv, iside, iclosed, gate, ...
            geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
            Rbound+rnear, Me, ones(K+1,1), Pe, Pe, cat(3,R{:}), [C{:}], 0);
  % naive SLP potential at Pe + ONE handle-based eval correction (source K*Nnod -> target Me),
  % replacing the whole per-particle closesize/closeslp/closecorrapply loop.
  u=Lap3dSLPfmm(struct('x',Pe), struct('x',Xall,'w',wall), sigma, fmmeps);
  u=axpso_corr_apply_mex(heval, K*Nnod, sigma, Me, u);       % arbitrary targets: nu=Me /= nx
  tevaloff=toc(teo);
  Ue=uexf(Pe).'; err=abs(u-Ue)/max(abs(Ue));

  % ---- in-loop visualization: field error at the check points, updated per K ----
  dsurf=zeros(1,Me);
  for i=1:Me, dsurf(i)=min(vecnorm(Xall-Pe(:,i))); end
  figure(2),clf; hold on;
  for k=1:K, plot3(X{k}(1,:),X{k}(2,:),X{k}(3,:),'.','Color',[.7 .7 .7],'MarkerSize',2); end
  scatter3(Pe(1,:),Pe(2,:),Pe(3,:),14,log10(max(err(:).',1e-17)),'filled');
  plot3(Y(1,:),Y(2,:),Y(3,:),'p','MarkerSize',10,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
  cb=colorbar; cb.Label.String='log_{10}|u_h-u_{exact}|/max|u_{exact}|';
  axis equal; view(35,18); grid on; colormap('jet');
  title(sprintf('K=%d field error at check points',K)); drawnow;
  % exportgraphics(figure(2),sprintf('axissymslap_lap_slpn_max_physmat_K%d.png',K),'Resolution',200)
  fprintf('  err by shell:  d~0.05: %.2e   d~0.3: %.2e   d~1.0: %.2e\n', ...
          max(err(dsurf<0.15)), max(err(dsurf>=0.15 & dsurf<0.6)), max(err(dsurf>=0.6)));

  % 8. the data point (two lines: SETUP = OpenMP solver-module pass (mex built with
  %    OMP=ON -- default libomp thread count = ncores); EVAL/SOLVE = FMM-dominated)
  dof=K*Nnod;                                                 % SCALAR: 1 dof per node
  fprintf(['K=%4d  dof=%8d :  T_setup %6.1fs (T_self %6.1fs + T_cross %6.1fs)', ...
           '   [%d threads]   %8.0f dof/s/core\n'], ...
          K, dof, tself+tcross, tself, tcross, ncores, dof/(tself+tcross)/ncores);
  fprintf(['                        T_eval %6.3f s/iter (fmm eps %.0e, %4.0f us/src)  iters %3d (flag=%d, relres=%.1e)  ', ...
           'T_solve %5.0fs  T_eval_off %5.1fs (%d tgts)  err %.2e   [%d cores]  %8.0f dof/s/core\n'], ...
          Teval, fmmeps, Teval/(K*Nnod)*1e6, iter(2), flag, relres, tsolve, tevaloff, Me, max(err), ...
          ncores, dof/Teval/ncores);
end



function y=applytimed(x, fmmeps, Xall, Nall, wall, hself, hcross)
% timed K-particle matrix-free (-1/2 I + S') matvec: global lfmm3d naive normal-derivative
% (self i~=j excluded) + the handle-based near correction (self then cross, each ACCUMULATING).
% SCALAR density -> NO lab<->particle R sandwich.  T_eval = accumulated apply time / count.
global APPLY_T APPLY_N
t0=tic;
srcinfo.sources=Xall; srcinfo.charges=(x(:).').*wall;        % scalar charges, weights baked in
U=lfmm3d(fmmeps,srcinfo,2);                                   % pot + grad AT sources (self i~=j excluded)
y=(sum(U.grad.*Nall,1)).';                                   % S' = grad . node normal, scalar per node
n=numel(x);
y=axpso_corr_apply_mex(hself,  n, x, n, y);                  % self  correction (accumulate; nu==nx)
y=axpso_corr_apply_mex(hcross, n, x, n, y);                  % cross correction (accumulate; nu==nx)
APPLY_T=APPLY_T+toc(t0); APPLY_N=APPLY_N+1;
end
function u=lapsum(P,Y,Q)
u=zeros(1,size(P,2));
for j=1:size(Y,2), u=u+Q(j)./(4*pi*vecnorm(P-Y(:,j))); end
end
