clearvars; format short e;
addpath('../../../axibie/utils');
addpath('../../utils');
addpath('../../matlab');
addpath('../../external/fmm3d/matlab');                       % stfmm3d (via Sto3d{SLP,DLP}fmm_il)

% MULTI-particle Stokes DLP exterior-Dirichlet BVP -- PHYSICAL-SPACE SPARSE-CORRECTION version of
% test_axissymsstok_stok_dlp_bvp_multi.m, in the framework style of the SLPn twin
% test_axissymsstok_stok_slpn_bvp_multi_physmat.m and the Laplace DLP twin
% test_axissymslap_lap_dlp_bvp_multi_physmat.m.  COMBINED-FIELD (D+S)[sigma]: the +1/2 I exterior jump
% rides in the DLP self block, the SLP regularizes the DLP rigid-body null space -> full-rank (no
% deflation).  Fully MATRIX-FREE mex; each near correction is FULL close-eval minus naive, in Fortran:
%   SETUP (per layer, per pass): axpso_close_setup_mex returns the FULL close-eval matrix (iform=0)
%     via the WHOLE-PARTICLE closestokdlp/closestokslp build, then axpso_close2corr_set_mex subtracts
%     the naive Sto3d kernel IN FORTRAN (compact layout, self node auto-zero) -> the handle holds the
%     correction.  Four SYSTEM handles: {D,S} x {self,cross}.
%   SYSTEM (velocity): Sto3dSLPfmm_il + Sto3dDLPfmm_il naive (D+S) at the sources (targets==sources ->
%     self dropped; the +1/2 jump rides in the DLP self correction) + four corr applies (ACCUMULATE,
%     the lab<->local R sandwich INSIDE the module).
%   SOLVE: restarted gmres, combined-field second-kind well-conditioned (no deflation).
%   EVAL (velocity): naive (D+S) FMM at the grid + {D,S} eval corrections, same recipe.
% VECTOR (nc=3), node-interleaved lab frame; mu=1.

p=16; np_vals=2:8; errmax=nan(1,numel(np_vals)); mu=1.0; fmmeps=1e-12; rfac=1.25;   % rfac: near-zone enlargement (stresslet far-Nv annulus)
rot=@(u,th) cos(th)*eye(3)+sin(th)*[0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0]+(1-cos(th))*(u*u');
Zf={@(t) 0.5*sin(t)-1i*1.0*cos(t), @(t) 1.0*sin(t)-1i*0.5*cos(t)}; ac=[0.5 1.0; 1.0 0.5];
stk=@(P,y,Fv)(1/(8*pi*mu))*( Fv./vecnorm(P-y) + (Fv.'*(P-y)).*(P-y)./vecnorm(P-y).^3 );    % interior Stokeslet velocity
yloc={[0.10 0.05; 0 0.05; 0.20 -0.30],[0.20 -0.10; 0 0.15; 0.10 -0.10]};
Floc={[1.0 -0.5; -0.7 0.6; 0.5 -0.4],[0.6 -0.4; 0.3 -0.5; -0.7 0.5]};

for kk=1:numel(np_vals)
  np=np_vals(kk); M=2*np; nphi=2*M+1; ph=2*pi*(0:nphi-1)/nphi; iside=1; iclosed=0;
  chi_crit=cosh(log(1e13)/M);                                  % resolution-based near-zone (== dense multi)
  rng(1);

  % 1. particles: random pose; LOCAL-frame meridian geometry + 3D ring nodes/normals (lab)
  C=cell(2,1); R=cell(2,1); X=cell(2,1); NX=cell(2,1); sq=cell(2,1); geo=cell(2,1);
  for k=1:2
    s=[]; s.p=p; s.Z=Zf{k}; s.tpan=linspace(0,pi,np+1)'; s=quadr(s,[],'p','G');
    d=randn(3,1); C{k}=(k-1)*1.6*d/norm(d); uu=randn(3,1); uu=uu/norm(uu); R{k}=rot(uu,2*pi*rand);
    rho=real(s.x); z=imag(s.x); nr=real(s.nx); nz=imag(s.nx);
    X{k}=R{k}*[(kron(cos(ph(:)),rho)).';(kron(sin(ph(:)),rho)).';(kron(ones(nphi,1),z)).']+C{k};
    NX{k}=R{k}*[(kron(cos(ph(:)),nr)).';(kron(sin(ph(:)),nr)).';(kron(ones(nphi,1),nz)).'];   % 3D node normals (lab)
    sq{k}=axisym_to_3d_quadrature(real(s.x),imag(s.x),s.ws,nphi);
    g=[]; g.sx=s.x(:); g.snx=s.nx(:); g.sws=s.ws(:); g.swxp=s.wxp(:); g.tpan=s.tpan(:); geo{k}=g;
  end
  N=numel(geo{1}.sx); Nnod=N*nphi; Ndof=3*Nnod;
  Xall=[X{1} X{2}]; Nall=[NX{1} NX{2}]; wall=[sq{1}.w sq{2}.w];

  % interior Stokeslets -> exact exterior velocity uex (Dirichlet data on Xall, node-interleaved)
  Y=zeros(3,0); Fg=zeros(3,0); for k=1:2, Y=[Y, R{k}*yloc{k}+C{k}]; Fg=[Fg, R{k}*Floc{k}]; end
  uex=@(P) stk(P,Y(:,1),Fg(:,1))+stk(P,Y(:,2),Fg(:,2))+stk(P,Y(:,3),Fg(:,3))+stk(P,Y(:,4),Fg(:,4));
  g_di=reshape(uex(Xall),[],1);                               % node-interleaved velocity data

  % 2. flat per-particle meridian geometry consumed by the handle-based setups (self + cross)
  K2=2;
  sxs=zeros(N,K2); snxs=sxs; swss=zeros(N,K2); swxps=zeros(N,K2); tpans=zeros(np+1,K2);
  for k=1:K2, gs=geo{k}; sxs(:,k)=gs.sx; snxs(:,k)=gs.snx; swss(:,k)=gs.sws; swxps(:,k)=gs.swxp; tpans(:,k)=gs.tpan; end
  pv=p*ones(K2,1); npv=np*ones(K2,1); pmv=M*ones(K2,1);
  geomoff=1+(0:K2)'*N; tpanoff=1+(0:K2)'*(np+1); targoff=1+(0:K2)'*Nnod;
  nsx=K2*N; ntpan=K2*(np+1);
  sxf=sxs(:); snxf=snxs(:); swsf=swss(:); swxpf=swxps(:); tpanf=tpans(:);
  parcs=arrayfun(@(j) sum(geo{1}.sws((j-1)*p+(1:p))), 1:np);       % panel arclengths
  rnear=rfac*max(max(real(geo{1}.sx))*sqrt(2*(chi_crit-1)), 2*max(parcs));  % rfac x dense multi's gs.rnear
  gate=rnear/sqrt(min(parcs));                                     % EVERY panel ball gate*sqrt(arc_j) >= rnear
  rball=1.0+gate*sqrt(max(parcs));                                 % canonical-ball bound >= largest panel ball

  % 3. SYSTEM corrections: FULL close-eval minus naive for BOTH layers, self and cross (D = stresslet
  %    with source normals, S = Stokeslet velocity).  hDself/hSself carry the +1/2 jump (DLP self).
  hDself = axpso_close_setup_mex(2, 3, 1, 0, mu, 1, K2, pv, npv, pmv, iside, iclosed, gate, ...
             geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
             rball, K2*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);         % DLP self
  axpso_close2corr_set_mex(hDself, K2, geomoff, nsx, sxf, snxf, swsf, targoff, K2*Nnod, Xall, Nall, [C{:}]);
  hSself = axpso_close_setup_mex(2, 1, 1, 0, mu, 1, K2, pv, npv, pmv, iside, iclosed, gate, ...
             geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
             rball, K2*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);         % SLP self
  axpso_close2corr_set_mex(hSself, K2, geomoff, nsx, sxf, snxf, swsf, targoff, K2*Nnod, Xall, Nall, [C{:}]);
  hDcross = axpso_close_setup_mex(2, 3, 1, 0, mu, 2, K2, pv, npv, pmv, iside, iclosed, gate, ...
              geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
              rball, K2*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);        % DLP cross
  axpso_close2corr_set_mex(hDcross, K2, geomoff, nsx, sxf, snxf, swsf, targoff, K2*Nnod, Xall, Nall, [C{:}]);
  hScross = axpso_close_setup_mex(2, 1, 1, 0, mu, 2, K2, pv, npv, pmv, iside, iclosed, gate, ...
              geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
              rball, K2*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);        % SLP cross
  axpso_close2corr_set_mex(hScross, K2, geomoff, nsx, sxf, snxf, swsf, targoff, K2*Nnod, Xall, Nall, [C{:}]);
  % peek: the DLP self correction block (full - naive)
  for k=1:K2
    gs=geo{k}; t3d0=[real(gs.sx).'; zeros(1,N); imag(gs.sx).'];
    tcxi0=zeros(np+1,1); nt=0;
    [~,nt]=axps_closesize_mex(N, gs.sx, t3d0, p, np, gs.sx, gs.sws, gate, tcxi0, nt);
    fprintf('  DLP self block k=%d:  %d x %d  (ntcx=%d)\n', k, 3*nt, 3*nphi*p, nt);
  end

  % 4. matrix-free operator: FMM naive (D+S) + four corr applies (R sandwich INSIDE the module)
  applyA=@(x) applymvDS(x, fmmeps, Xall, Nall, wall, hDself, hDcross, hSself, hScross);

  % 5. solve (restarted gmres, combined-field second-kind well-conditioned -> no deflation)
  [sigma,flag,relres,iter]=gmres(applyA, g_di, 50, 1e-11, 100);
  fprintf('  gmres: flag=%d  iters=%d  relres=%.3e\n', flag, iter(end), relres);

  % 6. 3D exterior target grid (outside both particles)
  cmid=(C{1}+C{2})/2; gv=linspace(-3,3,22); [GX,GY,GZ]=ndgrid(gv,gv,gv); P=[GX(:)';GY(:)';GZ(:)']+cmid;
  ext=true(1,size(P,2)); dmin=inf(1,size(P,2));
  for k=1:2, loc=R{k}.'*(P-C{k}); sd=(hypot(loc(1,:),loc(2,:))/ac(k,1)).^2+(loc(3,:)/ac(k,2)).^2; ext=ext&(sd>1.05); dmin=min(dmin,sd-1); end
  Pe=P(:,ext); Me=size(Pe,2); dmin=dmin(ext);

  % 7. field eval u = (D+S)[sigma]: FMM naive (D+S) at Pe + {D,S} eval corrections (targoff=ones ->
  %    every source sees ALL Me grid targets), same full-minus-naive recipe applied source 3*K*Nnod -> 3*Me.
  u = Sto3dSLPfmm_il(struct('x',Pe), struct('x',Xall,'w',wall), sigma, fmmeps) ...
    + Sto3dDLPfmm_il(struct('x',Pe), struct('x',Xall,'w',wall,'nx',Nall), sigma, fmmeps);
  hDeval = axpso_close_setup_mex(2, 3, 1, 0, mu, 2, K2, pv, npv, pmv, iside, iclosed, gate, ...
             geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
             rball, Me, ones(K2+1,1), Pe, Pe, cat(3,R{:}), [C{:}], 0);            % DLP eval
  axpso_close2corr_set_mex(hDeval, K2, geomoff, nsx, sxf, snxf, swsf, ones(K2+1,1), Me, Pe, Pe, [C{:}]);
  hSeval = axpso_close_setup_mex(2, 1, 1, 0, mu, 2, K2, pv, npv, pmv, iside, iclosed, gate, ...
             geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
             rball, Me, ones(K2+1,1), Pe, Pe, cat(3,R{:}), [C{:}], 0);            % SLP eval
  axpso_close2corr_set_mex(hSeval, K2, geomoff, nsx, sxf, snxf, swsf, ones(K2+1,1), Me, Pe, Pe, [C{:}]);
  u=axpso_corr_apply_mex(hDeval, 3*K2*Nnod, sigma, 3*Me, u);
  u=axpso_corr_apply_mex(hSeval, 3*K2*Nnod, sigma, 3*Me, u);
  U3=reshape(u,3,Me); err=vecnorm(U3-uex(Pe)).'/max(vecnorm(uex(Pe))); inb=dmin.'<0.1; errmax(kk)=max(err);
  fprintf('np=%2d  N=%5d/part  M=%2d  M3=%5d :  near(d<0.1) %.2e   far(d>=0.1) %.2e\n', ...
          np, Ndof, M, Me, max([err(inb);0]), max(err(~inb)));

  figure(2),clf; hold on;
  scatter3(Pe(1,:),Pe(2,:),Pe(3,:),14,log10(max(err,1e-17)),'filled');
  plot3(X{1}(1,:),X{1}(2,:),X{1}(3,:),'.','Color',[.6 .6 .6]);
  plot3(X{2}(1,:),X{2}(2,:),X{2}(3,:),'.','Color',[.6 .6 .6]);
  plot3(Y(1,:),Y(2,:),Y(3,:),'p','MarkerSize',16,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
  cb=colorbar; cb.Label.String='log_{10}|u_h-u_{exact}|/max|u_{exact}|';
  axis equal; view(35,18); grid on; colormap('jet');
  xlabel('x'); ylabel('y'); zlabel('z'); title('2-particle Stokes DLP (D+S, physmat sparse-correction): 3D velocity error (last refinement)');
  figure(2),ylim([-1 3])
end

figure(2),clf; hold on;
scatter3(Pe(1,:),Pe(2,:),Pe(3,:),14,log10(max(err,1e-17)),'filled');
plot3(X{1}(1,:),X{1}(2,:),X{1}(3,:),'.','Color',[.6 .6 .6]);
plot3(X{2}(1,:),X{2}(2,:),X{2}(3,:),'.','Color',[.6 .6 .6]);
plot3(Y(1,:),Y(2,:),Y(3,:),'p','MarkerSize',16,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
cb=colorbar; cb.Label.String='log_{10}|u_h-u_{exact}|/max|u_{exact}|';
clim([-16 -4]); axis equal; view(35,18); grid on; colormap('jet');
xlabel('x'); ylabel('y'); zlabel('z'); title('2-particle Stokes DLP (D+S, physmat sparse-correction): 3D velocity error (last refinement)');
figure(2),ylim([-1 3])

% --- convergence plot: max field error vs panel count ---
figure(1),clf; semilogy(np_vals,errmax,'o-k'); grid on; xlabel('n_p'); ylabel('max field err');
title('multi-particle Stokes DLP (D+S, physmat sparse-correction): h-refinement, p=16, M=2 n_p');

% cmp = getPyPlot_cMap('rainbow', [], [], '"/Users/hzhu/.pyenv/versions/3.11.13/bin/python"');
% colormap(cmp)

% exportgraphics(figure(1),'axissymsstok_stok_dlp_multi_physmat_convergence.png','Resolution',200)
% exportgraphics(figure(2),'axissymsstok_stok_dlp_multi_physmat_error.png','Resolution',200)

function y = applymvDS(x, fmmeps, Xall, Nall, wall, hDself, hDcross, hSself, hScross)
% matrix-free (+1/2 I + D + S) matvec: naive (D+S) = Sto3dSLPfmm_il + Sto3dDLPfmm_il at the sources
% (targets==sources==Xall -> self dropped; +1/2 jump in the DLP self correction) + four handle corr
% applies (ACCUMULATE, R sandwich in-module).  x node-interleaved (3 per node).
y = Sto3dSLPfmm_il(struct('x',Xall), struct('x',Xall,'w',wall), x, fmmeps) ...
  + Sto3dDLPfmm_il(struct('x',Xall), struct('x',Xall,'w',wall,'nx',Nall), x, fmmeps);
n=numel(x);
y=axpso_corr_apply_mex(hDself,  n, x, n, y);
y=axpso_corr_apply_mex(hDcross, n, x, n, y);
y=axpso_corr_apply_mex(hSself,  n, x, n, y);
y=axpso_corr_apply_mex(hScross, n, x, n, y);
end
