clearvars; format short e;
addpath('../../../axibie/utils');
addpath('../../utils');
addpath('../../matlab');
addpath('../../external/fmm3d/matlab');                       % stfmm3d

% K-particle Stokes DLP combined-field (D+S) exterior-Dirichlet BVP -- the Dirichlet sibling of
% test_axissymsstok_stok_slpn_bvp_max_physmat5.m (same lattice / timing / grid-eval conventions),
% generalizing test_axissymsstok_stok_dlp_bvp_multi_physmat.m from 2 to K particles.
% COMBINED-FIELD (D+S)[sigma] = u_ex: the +1/2 I jump rides in the DLP self correction, the SLP
% regularizes the DLP rigid-body null space -> SECOND-KIND full-rank; gmres needs NO deflation and
% NO preconditioner (iteration count ~flat in K -- the ideal scaling operator).
%   SETUP: FOUR corr_setup handles {D,S} x {self,cross}; the D and S corrections share the SAME
%     topology (same gate/targets), so they are COMBINED per particle -- corr_get(D) + corr_get(S)
%     -> corr_set into the D handle -- and the matvec applies ONE combined handle per pass (halves
%     the apply cost; the donated S handles remain allocated -- the framework has no release verb).
%     ntcx per particle is recovered by axps_closesize over a candidate SUPERSET (exact count: far
%     candidates fail the per-panel gate).
%   SYSTEM (velocity): ONE stfmm3d call (stoklet + strslet/strsvec, pot at sources, self i~=j
%     excluded) = naive (D+S) + the two combined corr applies (R sandwich INSIDE the module).
%   EVAL: stfmm3d (D+S) pottarg at the grid + ONE combined {D,S} eval handle apply.
% Near-zone: gate=2.0 (the max-family convention): at the fixed np=4 the discretization floor
% (~1e-6) sits far above the stresslet far-Nv annulus, so the multi test's rfac/chi_crit deep-
% refinement knob is not needed here.

global APPLY_T APPLY_N                                        % T_eval accumulator (set in applytimed)

p=16; np=4; mu=1.0; gate=2.0; iside=1; iclosed=0;
M=2*np; nphi=2*M+1;
Ksides=[2 3 4];                                               % K = 8, 27, 64
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
ncores=feature('numcores');
fprintf('p=%d np=%d nphi=%d (Nnod=%d/particle)  gate=%.1f  rnear=%.2f  lattice spacing %.2f  ncores=%d\n', ...
        p, np, nphi, Nnod, gate, rnear, Lsp, ncores);

% self ntcx per SHAPE (targets = azimuth-0 meridian reps; superset-exact, shared by all particles
% of that shape) -- needed to corr_get/corr_set the combine below
ntS=zeros(2,1);
for sh=1:2
  gs=geo{sh}; t3d0=[real(gs.sx).'; zeros(1,N); imag(gs.sx).'];
  tcxi0=zeros(np+1,1); nt=0;
  [~,nt]=axps_closesize_mex(N, gs.sx, t3d0, p, np, gs.sx, gs.sws, gate, tcxi0, nt);
  ntS(sh)=nt;
end

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
  for k=1:K
    nodeidx=(k-1)*Nnod+(1:Nnod);
    Xall(:,nodeidx)=X{k}; Nall(:,nodeidx)=NX{k}; wall(nodeidx)=sq{ty(k)}.w;
  end

  % 2. interior Stokeslets (2 per particle) -> exact exterior velocity (Dirichlet data)
  Y=zeros(3,2*K); Fg=zeros(3,2*K);
  for k=1:K, Y(:,2*k-1:2*k)=R{k}*yloc{ty(k)}+C{k}; Fg(:,2*k-1:2*k)=R{k}*Floc{ty(k)}; end
  uexf=@(P) stksum(P,Y,Fg,stk);
  g_di=uexf(Xall); g_di=g_di(:);                              % exact velocity, node-interleaved (lab)

  % flat per-particle meridian geometry (consumed by all corr_setup calls)
  sxs=zeros(N,K); snxs=sxs; swss=zeros(N,K); swxps=zeros(N,K); tpans=zeros(np+1,K);
  for k=1:K
    gs=geo{ty(k)};
    sxs(:,k)=gs.sx; snxs(:,k)=gs.snx; swss(:,k)=gs.sws; swxps(:,k)=gs.swxp; tpans(:,k)=gs.tpan;
  end
  pv=p*ones(K,1); npv=np*ones(K,1); pmv=M*ones(K,1);
  geomoff=1+(0:K)'*N; tpanoff=1+(0:K)'*(np+1); targoff=1+(0:K)'*Nnod;
  nsx=K*N; ntpan=K*(np+1);
  sxf=sxs(:); snxf=snxs(:); swsf=swss(:); swxpf=swxps(:); tpanf=tpans(:);

  % 3. SELF corrections {D,S} + COMBINE: two corr_setup calls (ilayer=3 DLP stresslet with the
  %    +1/2 jump, ilayer=1 SLP Stokeslet), then per particle corr_get(D)+corr_get(S) ->
  %    corr_set into the D handle; the matvec applies ONE combined handle.
  tself=tic;
  [hDself, nbDs] = axpso_corr_setup_mex(2, 3, 1, 0, mu, 1, K, pv, npv, pmv, iside, iclosed, gate, ...
            geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
            Rbound+rnear, K*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);
  [hSself, nbSs] = axpso_corr_setup_mex(2, 1, 1, 0, mu, 1, K, pv, npv, pmv, iside, iclosed, gate, ...
            geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
            Rbound+rnear, K*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);
  valsDS=[]; iksDS=[];
  for k=1:K
    nt=ntS(ty(k));
    [bD,ikk,~]=axpso_corr_get_mex(hDself,k,nt,3,p,np,M,zeros(9*nt*nphi*p,1),zeros(nt,1),zeros(np+1,1));
    [bS,~,~]=axpso_corr_get_mex(hSself,k,nt,3,p,np,M,zeros(9*nt*nphi*p,1),zeros(nt,1),zeros(np+1,1));
    valsDS=[valsDS; bD+bS]; iksDS=[iksDS; ikk];
  end
  axpso_corr_set_mex(hDself, valsDS, iksDS);    % hDself now holds the (D+S) correction
  tself=toc(tself); selfMB=(nbDs+nbSs)/1e6;

  % 4. CROSS corrections {D,S} + COMBINE: same recipe (iinter=2); cross ntcx per particle is
  %    recovered by closesize over ALL non-own nodes rotated into the source frame (superset-exact).
  tcross=tic;
  [hDcross, nbDx] = axpso_corr_setup_mex(2, 3, 1, 0, mu, 2, K, pv, npv, pmv, iside, iclosed, gate, ...
             geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
             Rbound+rnear, K*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);
  [hScross, nbSx] = axpso_corr_setup_mex(2, 1, 1, 0, mu, 2, K, pv, npv, pmv, iside, iclosed, gate, ...
             geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
             Rbound+rnear, K*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);
  valsDS=[]; iksDS=[];
  for k=1:K
    gs=geo{ty(k)};
    other=setdiff(1:K*Nnod,(k-1)*Nnod+(1:Nnod)); loc=R{k}.'*(Xall(:,other)-C{k});
    ttx=complex(hypot(loc(1,:),loc(2,:)).',loc(3,:).'); tc2=zeros(np+1,1); nt=0;
    [~,nt]=axps_closesize_mex(numel(other), ttx, loc, p, np, gs.sx, gs.sws, gate, tc2, nt);
    if nt==0, continue; end
    [bD,ikk,~]=axpso_corr_get_mex(hDcross,k,nt,3,p,np,M,zeros(9*nt*nphi*p,1),zeros(nt,1),zeros(np+1,1));
    [bS,~,~]=axpso_corr_get_mex(hScross,k,nt,3,p,np,M,zeros(9*nt*nphi*p,1),zeros(nt,1),zeros(np+1,1));
    valsDS=[valsDS; bD+bS]; iksDS=[iksDS; ikk];
  end
  axpso_corr_set_mex(hDcross, valsDS, iksDS);   % hDcross now holds the (D+S) correction
  tcross=toc(tcross); crossMB=(nbDx+nbSx)/1e6;

  % 5. matrix-free operator (timed per apply): ONE stfmm3d (D+S) + TWO combined corr applies.
  %    SECOND-KIND full-rank -> no deflation, no preconditioner.
  APPLY_T=0; APPLY_N=0;
  applyA=@(x) applytimed(x, fmmeps, Xall, Nall, wall, hDself, hDcross);

  % 6. solve;  T_eval = mean apply time over ALL applies inside the solve
  tso=tic;
  [sigma,flag,relres,iter]=gmres(applyA, g_di, [], 1e-9, 200);
  tsolve=toc(tso);
  Teval=APPLY_T/APPLY_N;

  % 7. off-surface eval: regular ng^3 grid, exact spheroid inside-filter; stfmm3d (D+S)
  %    velocity at the outside targets + ONE combined {D,S} eval handle apply.
  teo=tic;
  ng=16; margin=0.9;
  glo=-Lsp*(Kside-1)/2 - Rbound - margin;
  xs=linspace(glo,-glo,ng);
  [xg,yg,zg]=ndgrid(xs,xs,xs);
  Pg=[xg(:).'; yg(:).'; zg(:).'];
  inside=false(1,size(Pg,2));
  for k=1:K
    loc=R{k}.'*(Pg-C{k});
    sd=(hypot(loc(1,:),loc(2,:))/ac(ty(k),1)).^2+(loc(3,:)/ac(ty(k),2)).^2;
    inside=inside | (sd<1.0);
  end
  Pe=Pg(:,~inside); Me=size(Pe,2);
  fprintf('  targets: %d^3 grid, %d outside particles (h=%.3f)\n', ng, Me, xs(2)-xs(1));
  sig3=reshape(sigma,3,[]).*wall;
  srcinfo=[]; srcinfo.sources=Xall; srcinfo.stoklet=sig3; srcinfo.strslet=sig3; srcinfo.strsvec=Nall;
  U=stfmm3d(fmmeps,srcinfo,0,Pe,1);
  uc=U.pottarg(:);                                            % naive (D+S) velocity, node-interleaved
  hDeval = axpso_corr_setup_mex(2, 3, 1, 0, mu, 2, K, pv, npv, pmv, iside, iclosed, gate, ...
            geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
            Rbound+rnear, Me, ones(K+1,1), Pe, Pe, cat(3,R{:}), [C{:}], 0);
  hSeval = axpso_corr_setup_mex(2, 1, 1, 0, mu, 2, K, pv, npv, pmv, iside, iclosed, gate, ...
            geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
            Rbound+rnear, Me, ones(K+1,1), Pe, Pe, cat(3,R{:}), [C{:}], 0);
  valsDS=[]; iksDS=[];
  for k=1:K                                                   % combine {D,S} eval (targets = ALL Pe, superset-exact)
    gs=geo{ty(k)};
    loc=R{k}.'*(Pe-C{k});
    ttx=complex(hypot(loc(1,:),loc(2,:)).',loc(3,:).'); tc2=zeros(np+1,1); nt=0;
    [~,nt]=axps_closesize_mex(Me, ttx, loc, p, np, gs.sx, gs.sws, gate, tc2, nt);
    if nt==0, continue; end
    [bD,ikk,~]=axpso_corr_get_mex(hDeval,k,nt,3,p,np,M,zeros(9*nt*nphi*p,1),zeros(nt,1),zeros(np+1,1));
    [bS,~,~]=axpso_corr_get_mex(hSeval,k,nt,3,p,np,M,zeros(9*nt*nphi*p,1),zeros(nt,1),zeros(np+1,1));
    valsDS=[valsDS; bD+bS]; iksDS=[iksDS; ikk];
  end
  axpso_corr_set_mex(hDeval, valsDS, iksDS);    % (D+S) via ONE-call whole-record set
  uc=axpso_corr_apply_mex(hDeval, 3*K*Nnod, sigma, 3*Me, uc);
  Uh=zeros(3,size(Pg,2)); Uh(:,~inside)=reshape(uc,3,[]);
  tevaloff=toc(teo);
  Ue=uexf(Pe);
  scal=max(vecnorm(Ue));
  errg=zeros(1,size(Pg,2)); errg(~inside)=vecnorm(Uh(:,~inside)-Ue)/scal;
  err=max(errg);

  % ---- in-loop visualization: field error at the check points, updated per K ----
  errp=errg(~inside);                                         % per-point error at the Pe targets
  dsurf=zeros(1,Me);
  for i=1:Me, dsurf(i)=min(vecnorm(Xall-Pe(:,i))); end        % dist to nearest surface node
  figure(2),clf; hold on;
  for k=1:K, plot3(X{k}(1,:),X{k}(2,:),X{k}(3,:),'.','Color',[.7 .7 .7],'MarkerSize',2); end
  scatter3(Pe(1,:),Pe(2,:),Pe(3,:),14,log10(max(errp,1e-17)),'filled');
  plot3(Y(1,:),Y(2,:),Y(3,:),'p','MarkerSize',10,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
  cb=colorbar; cb.Label.String='log_{10}|u_h-u_{exact}|/max|u_{exact}|';
  axis equal; view(35,18); grid on; colormap('jet');
  title(sprintf('K=%d field error at check points',K)); drawnow;
  % exportgraphics(figure(2),sprintf('axissymsstok_stok_slpn_max_physmat_K%d.png',K),'Resolution',200)
  fprintf('  err by shell:  d~0.05: %.2e   d~0.3: %.2e   d~1.0: %.2e\n', ...
          max(errp(dsurf<0.15)), max(errp(dsurf>=0.15 & dsurf<0.6)), max(errp(dsurf>=0.6)));

  % 8. the data point (paper-style, two lines, == physmat5 conventions)
  dof=3*K*Nnod;
  fprintf(['K=%4d  dof=%8d :  T_setup %6.1fs (T_self %6.1fs %.0f MB + T_cross %6.1fs %.0f MB, {D,S} combined)', ...
           '   [%d threads]   %8.0f dof/s/core\n'], ...
          K, dof, tself+tcross, tself, selfMB, tcross, crossMB, ncores, dof/(tself+tcross)/ncores);
  fprintf(['                        T_eval %6.3f s/iter (fmm eps %.0e, %4.0f us/src)  iters %3d (flag=%d, relres=%.1e)  ', ...
           'T_solve %5.0fs  T_eval_off %5.1fs (%d tgts)  err %.2e   [%d cores]  %8.0f dof/s/core\n'], ...
          Teval, fmmeps, Teval/(K*Nnod)*1e6, iter(2), flag, relres, tsolve, tevaloff, Me, max(err), ...
          ncores, dof/Teval/ncores);

end

function y=applytimed(x, fmmeps, Xall, Nall, wall, hDSself, hDScross)
% timed K-particle matrix-free (+1/2 I + D + S) matvec: ONE stfmm3d call (stoklet + strslet) =
% naive (D+S) velocity at the sources (self i~=j excluded; the +1/2 jump rides in the combined
% self correction) + the two COMBINED {D,S} handle applies (R sandwich INSIDE the module).
global APPLY_T APPLY_N
t0=tic;
sig3=reshape(x,3,[]).*wall;
srcinfo.sources=Xall; srcinfo.stoklet=sig3; srcinfo.strslet=sig3; srcinfo.strsvec=Nall;
U=stfmm3d(fmmeps,srcinfo,1);
y=U.pot(:);
n=numel(x);
y=axpso_corr_apply_mex(hDSself,  n, x, n, y);                 % combined (D+S) self  correction
y=axpso_corr_apply_mex(hDScross, n, x, n, y);                 % combined (D+S) cross correction
APPLY_T=APPLY_T+toc(t0); APPLY_N=APPLY_N+1;
end
function u=stksum(P,Y,Fg,stk)
u=zeros(3,size(P,2));
for j=1:size(Y,2), u=u+stk(P,Y(:,j),Fg(:,j)); end
end
