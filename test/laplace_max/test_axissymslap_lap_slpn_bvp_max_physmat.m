clearvars; format short e;
addpath('../../../axibie/utils');
addpath('../../utils');
addpath('../../matlab');
addpath('../../external/fmm3d/matlab');                       % lfmm3d / Lap3dSLPfmm
addpath('../../external/kdtree/toolbox');                    % KDTree (one canonical tree per K)

% K-particle Laplace SLPn (S') exterior-Neumann BVP -- the laplace_max SCALING BENCHMARK,
% the scalar analog of test_axissymsstok_stok_slpn_bvp_max_physmat.m (stokes_max, issue #43).
% K particles on a regular Kside^3 lattice (random rotations, lattice spacing chosen so the
% worst-case surface gap is 0.6 -- the validated near-zone regime).  Fully matrix-free:
%   naive S' = lfmm3d normal-derivative at the sources (self excluded; -1/2 I rides in the
%   self corrections) + PER-PARTICLE self corrections + CROSS corrections per NEIGHBOR pair
%   (centroid prefilter), GMRES.  NO deflation: the Laplace exterior-Neumann operator
%   (-1/2 + S') is second-kind FULL-RANK (cf. test_axissymslap_lap_slpn_bvp_multi.m) --
%   unlike Stokes, whose traction operator carries the per-particle normal null vectors.
%
% SCALAR (nc=1) simplifications vs the Stokes template:
%   * density sigma is ONE value per node (not a 3-vector) -> no lab<->particle R sandwich
%     on the density (a scalar is rotation-invariant);
%   * no viscosity mu; the FMM matvec is lfmm3d grad . normal instead of the stresslet.
% Close-correction mexes (axps sparse module -- Laplace lives in the SAME module):
%   axps_closesize_mex     EXISTS (kernel-independent near-panel sizing)
%   axps_closeslp_mex      EXISTS (Laplace SLP value close pass; used by the field eval)
%   axps_closecorrapply_mex EXISTS (nc-generic apply; called with nc=1)
%   axps_closeslpn_mex     MISSING -- the Laplace S' close pass, the one routine this
%                          benchmark needs implemented (scalar analog of axps_closestokslpn_mex).
% Timing convention (Li-Malhotra-Veerapaneni): T_self/T_cross = per-configuration operator
% SETUP; T_eval = per-GMRES-iteration apply time; T_eval_off = off-surface field eval.

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
  kd=KDTree(Xall.');                                          % ONE canonical 3D target tree

  % 2. interior point charges (2 per particle) -> exact exterior potential + surface flux data
  Y=zeros(3,2*K); Q=zeros(1,2*K);
  for k=1:K, Y(:,2*k-1:2*k)=R{k}*yloc{ty(k)}+C{k}; Q(2*k-1:2*k)=qloc{ty(k)}; end
  uexf=@(P) lapsum(P,Y,Q);
  g_tr=zeros(1,K*Nnod);                                       % dn(uex) = -sum_j q_j (n.(X-Y_j))/(4pi r^3)
  for j=1:2*K, d=Xall-Y(:,j); r=vecnorm(d); g_tr=g_tr - Q(j)*(sum(Nall.*d,1))./(4*pi*r.^3); end
  g_tr=g_tr(:);                                               % exact surface flux, scalar per node

  % 3. SELF corrections PER PARTICLE (generic path: own pass each, azimuth-0 orbit reps)
  tself=tic; Ssl=cell(K,1); isl=cell(K,1); tcxis=cell(K,1); ntcxs=zeros(K,1);
  for k=1:K
    gs=geo{ty(k)};
    t3d0=[real(gs.sx).'; zeros(1,N); imag(gs.sx).']; tn3d0=[real(gs.snx).'; zeros(1,N); imag(gs.snx).'];
    tcxi=zeros(np+1,1); ntcx=0;
    [tcxi,ntcx]=axps_closesize_mex(N, gs.sx, t3d0, p, np, gs.sx, gs.sws, gate, tcxi, ntcx);
    Sk=zeros(ntcx,nphi*p); ik=zeros(ntcx,1);                  % SCALAR block sizes (nc=1)
    [Sk,ik]=axps_closeslpn_mex(N, gs.sx, t3d0, tn3d0, p, np, nphi, gs.sx, gs.snx, gs.sws, gs.swxp, gs.tpan, gate, ...
        sq{ty(k)}.x, sq{ty(k)}.nx, sq{ty(k)}.w, M, iside, iclosed, ntcx, tcxi, Sk, ik);
    Ssl{k}=Sk; isl{k}=ik; tcxis{k}=tcxi; ntcxs(k)=ntcx;
  end
  tself=toc(tself); selfMB=sum(cellfun(@numel,Ssl))*8/1e6;

  % 4. CROSS corrections: ONE loop over source particles; targets = canonical nodes minus own block
  tcross=tic; Sx=cell(K,1); tcxix=cell(K,1); ntcxx=zeros(K,1);
  idxc=cell(K,1); canonl=cell(K,1); nu=zeros(K,1);
  for k=1:K
    cand = kd.ball(C{k}.', Rbound+rnear); cand=cand(:).';
    cand = cand(cand<=(k-1)*Nnod | cand>k*Nnod);             % drop own block (self is the circulant pass)
    if isempty(cand), continue, end
    loc = R{k}.'*(Xall(:,cand)-C{k}); nloc = R{k}.'*Nall(:,cand);   % geometry into particle frame
    ttx = (hypot(loc(1,:),loc(2,:)) + 1i*loc(3,:)).';
    gs=geo{ty(k)};
    tcxi=zeros(np+1,1); ntcx=0;
    [tcxi,ntcx]=axps_closesize_mex(numel(cand), ttx, loc, p, np, gs.sx, gs.sws, gate, tcxi, ntcx);
    ntcxx(k)=ntcx; tcxix{k}=tcxi;
    if ntcx==0, continue, end
    Sk=zeros(ntcx,nphi*p); ik=zeros(ntcx,1);
    [Sk,ik]=axps_closeslpn_mex(numel(cand), ttx, loc, nloc, p, np, nphi, gs.sx, gs.snx, gs.sws, gs.swxp, gs.tpan, gate, ...
        sq{ty(k)}.x, sq{ty(k)}.nx, sq{ty(k)}.w, M, iside, iclosed, ntcx, tcxi, Sk, ik);
    Sx{k}=Sk;
    canon = cand(ik).';                                       % local -> canonical node ids
    [canonl{k},~,idxc{k}] = unique(canon); nu(k)=numel(canonl{k});
  end
  npair=nnz(ntcxx); tcross=toc(tcross);
  crossMB=sum(cellfun(@numel,Sx))*8/1e6;

  % 5. matrix-free operator (timed per apply).  NO deflation: (-1/2 + S') is second-kind
  %    FULL-RANK for the Laplace exterior-Neumann problem (unlike the Stokes traction
  %    operator, whose per-particle normal null vectors need the rank-K deflator).
  APPLY_T=0; APPLY_N=0;
  applyA=@(x) applytimed(x, fmmeps, Xall, Nall, wall, blk, N, p, np, nphi, Nnod, ...
                         Ssl, isl, tcxis, ntcxs, Sx, idxc, canonl, nu, tcxix, ntcxx);

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
  u=Lap3dSLPfmm(struct('x',Pe), struct('x',Xall,'w',wall), sigma, fmmeps);   % naive SLP potential
  for k=1:K                                                   % per-particle SLP-value close correction
    gs=geo{ty(k)}; loc=R{k}.'*(Pe-C{k});
    if min(vecnorm(loc)) > Rbound+rnear, continue, end
    ttx=(hypot(loc(1,:),loc(2,:)) + 1i*loc(3,:)).';
    tcxi=zeros(np+1,1); ntcx=0;
    [tcxi,ntcx]=axps_closesize_mex(Me, ttx, loc, p, np, gs.sx, gs.sws, gate, tcxi, ntcx);
    if ntcx==0, continue, end
    Sk=zeros(ntcx,nphi*p); ik=zeros(ntcx,1);
    [Sk,ik]=axps_closeslp_mex(Me, ttx, loc, p, np, nphi, gs.sx, gs.snx, gs.sws, gs.swxp, gs.tpan, gate, ...
        sq{ty(k)}.x, sq{ty(k)}.nx, sq{ty(k)}.w, M, iside, iclosed, ntcx, tcxi, Sk, ik);
    uc=zeros(Me,1);
    uc=axps_closecorrapply_mex(1, Me, p, np, nphi, ntcx, tcxi, ik, Sk, sigma(blk{k}), 0, uc);
    u=u+uc;
  end
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
  fprintf('  err by shell:  d~0.05: %.2e   d~0.3: %.2e   d~1.0: %.2e\n', ...
          max(err(dsurf<0.15)), max(err(dsurf>=0.15 & dsurf<0.6)), max(err(dsurf>=0.6)));

  % 8. the data point (two lines: SETUP = serial single-core mex loop; EVAL/SOLVE = FMM-dominated)
  dof=K*Nnod;                                                 % SCALAR: 1 dof per node
  fprintf(['K=%4d  dof=%8d :  T_setup %6.1fs (T_self %6.1fs %.0f MB + T_cross %6.1fs %4d srcs %.0f MB)', ...
           '   [1 core]   %8.0f dof/s/core\n'], ...
          K, dof, tself+tcross, tself, selfMB, tcross, npair, crossMB, dof/(tself+tcross));
  fprintf(['                        T_eval %6.3f s/iter (fmm eps %.0e, %4.0f us/src)  iters %3d (flag=%d, relres=%.1e)  ', ...
           'T_solve %5.0fs  T_eval_off %5.1fs (%d tgts)  err %.2e   [%d cores]  %8.0f dof/s/core\n'], ...
          Teval, fmmeps, Teval/(K*Nnod)*1e6, iter(2), flag, relres, tsolve, tevaloff, Me, max(err), ...
          ncores, dof/Teval/ncores);
end

function y=applytimed(x, fmmeps, varargin)
% timed operator apply: accumulates total apply time / count for T_eval
global APPLY_T APPLY_N
t0=tic;
y=mvslpnK(x, fmmeps, varargin{:});
APPLY_T=APPLY_T+toc(t0); APPLY_N=APPLY_N+1;
end
function u=lapsum(P,Y,Q)
u=zeros(1,size(P,2));
for j=1:size(Y,2), u=u+Q(j)./(4*pi*vecnorm(P-Y(:,j))); end
end
function t=mvslpnK(x, fmmeps, Xall, Nall, wall, blk, N, p, np, nphi, Nnod, ...
                   Ssl, isl, tcxis, ntcxs, Sx, idxc, canonl, nu, tcxix, ntcxx)
% K-particle matrix-free (-1/2 I + S') matvec: global lfmm3d naive normal-derivative (self
% i~=j excluded) + PER-PARTICLE circulant self corrections + per-neighbor-pair cross
% corrections.  SCALAR density -> NO lab<->particle R sandwich (a scalar is rotation-invariant).
srcinfo.sources=Xall; srcinfo.charges=(x(:).').*wall;        % scalar charges, weights baked in
U=lfmm3d(fmmeps,srcinfo,2);                                   % pot + grad AT sources (self i~=j excluded)
t=(sum(U.grad.*Nall,1)).';                                    % S' = grad . node normal, scalar per node
K=numel(blk);
for k=1:K                                                     % SELF (per particle)
  us=zeros(Nnod,1);
  us=axps_closecorrapply_mex(1, N, p, np, nphi, ntcxs(k), tcxis{k}, isl{k}, Ssl{k}, x(blk{k}), 1, us);
  t(blk{k})=t(blk{k})+us;
end
for k=1:K                                                     % CROSS (canonical targets, compact)
  if ntcxx(k)==0, continue, end
  uc=zeros(nu(k),1);
  uc=axps_closecorrapply_mex(1, nu(k), p, np, nphi, ntcxx(k), tcxix{k}, idxc{k}, Sx{k}, x(blk{k}), 0, uc);
  t(canonl{k})=t(canonl{k})+uc;                              % scatter to canonical node ids (no rotation)
end
end
