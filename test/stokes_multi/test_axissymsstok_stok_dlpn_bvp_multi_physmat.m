clearvars; format short e;
addpath('../../../axibie/utils');
addpath('../../utils');
addpath('../../matlab');
addpath('../../external/fmm3d/matlab');                       % stfmm3d (+ Sto3d{SLP,DLP}fmm_il)

% MULTI-particle Stokes DLPn exterior-Neumann BVP -- PHYSICAL-SPACE SPARSE-CORRECTION version of
% test_axissymsstok_stok_dlpn_bvp_multi.m, in the framework style of the Laplace DLPn physmat twin.
% COMBINED-FIELD u = (S+D)[sigma]; match the surface traction (S'+D')[sigma] = t_ex (the -1/2 jump
% rides in the S' self correction; the S' regularizes the hypersingular D' -> full-rank).  Fully
% MATRIX-FREE mex; each near correction is FULL close-eval minus naive, in Fortran:
%   SETUP (per layer, per pass): axpso_close_setup_mex returns the FULL close-eval matrix (iform=0);
%     axpso_close2corr_set_mex subtracts the naive Sto3d kernel IN FORTRAN.  Four SYSTEM handles:
%     {S',D'} x {self,cross}.
%   SYSTEM (traction): ONE stfmm3d call (stoklet + strslet/strsvec) -> traction of (S+D) at the
%     sources (self i~=j excluded) + the four corr applies (ACCUMULATE, R sandwich in-module).
%   SOLVE: CALDERON right-preconditioned gmres -- (-1/2+S'+D') is NOT second-kind (||D'|| ~ 1/h^2,
%     sol err ~ kappa*relres unpreconditioned); right-precondition by the RANK-COMPLETED on-surface
%     S value apply P = S + sum_k nh_k nh_k' (bare S has the pressure null space S[n_k]=0, so
%     sigma = S tau could not carry the density's normal components); D'S = K'^2 - 1/4 I -> the
%     composite is second-kind, iterations flat.  sigma = P tau.
%   EVAL (velocity): Sto3dSLPfmm_il + Sto3dDLPfmm_il naive (S+D) at the grid + {S,D} eval corrections.
% VECTOR (nc=3), node-interleaved lab frame; mu=1.  Near-zone: chi_crit-based rnear x rfac with every
% panel ball >= rnear (the stresslet-family far-Nv lesson from the DLP physmat twin).

p=16; np_vals=2:8; errmax=nan(1,numel(np_vals)); mu=1.0; fmmeps=1e-12; rfac=1.25;
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
  N=numel(geo{1}.sx); Nnod=N*nphi; Ndof=3*Nnod; Mtot=2*Nnod;
  Xall=[X{1} X{2}]; Nall=[NX{1} NX{2}]; wall=[sq{1}.w sq{2}.w];
  blk={1:Ndof, Ndof+1:2*Ndof};

  % interior Stokeslets -> exact exterior velocity uex (eval) AND exact surface traction (Neumann data)
  Y=zeros(3,0); Fg=zeros(3,0); for k=1:2, Y=[Y, R{k}*yloc{k}+C{k}]; Fg=[Fg, R{k}*Floc{k}]; end
  uex=@(P) stk(P,Y(:,1),Fg(:,1))+stk(P,Y(:,2),Fg(:,2))+stk(P,Y(:,3),Fg(:,3))+stk(P,Y(:,4),Fg(:,4));
  [~,SpY]=Sto3dSLPmat_il(struct('x',Xall,'nx',Nall), struct('x',Y,'w',ones(1,size(Y,2))));
  g_tr = SpY*Fg(:);                                            % exact surface traction, node-interleaved (lab)

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

  % 3. SYSTEM corrections: FULL close-eval minus naive for BOTH traction layers, self and cross.
  %    The -1/2 jump rides in the S' self correction; D' = closestokdlpn (whole-particle split-open).
  hSpself = axpso_close_setup_mex(2, 2, 1, 0, mu, 1, K2, pv, npv, pmv, iside, iclosed, gate, ...
              geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
              rball, K2*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);            % SLPn self
  axpso_close2corr_set_mex(hSpself, K2, geomoff, nsx, sxf, snxf, swsf, targoff, K2*Nnod, Xall, Nall, [C{:}]);
  hDpself = axpso_close_setup_mex(2, 4, 1, 0, mu, 1, K2, pv, npv, pmv, iside, iclosed, gate, ...
              geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
              rball, K2*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);            % DLPn self
  axpso_close2corr_set_mex(hDpself, K2, geomoff, nsx, sxf, snxf, swsf, targoff, K2*Nnod, Xall, Nall, [C{:}]);
  hSpcross = axpso_close_setup_mex(2, 2, 1, 0, mu, 2, K2, pv, npv, pmv, iside, iclosed, gate, ...
               geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
               rball, K2*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);           % SLPn cross
  axpso_close2corr_set_mex(hSpcross, K2, geomoff, nsx, sxf, snxf, swsf, targoff, K2*Nnod, Xall, Nall, [C{:}]);
  hDpcross = axpso_close_setup_mex(2, 4, 1, 0, mu, 2, K2, pv, npv, pmv, iside, iclosed, gate, ...
               geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
               rball, K2*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);           % DLPn cross
  axpso_close2corr_set_mex(hDpcross, K2, geomoff, nsx, sxf, snxf, swsf, targoff, K2*Nnod, Xall, Nall, [C{:}]);

  % 4. matrix-free operator: stfmm3d naive (S'+D') traction + four corr applies
  applyA=@(x) applymvSD(x, fmmeps, Xall, Nall, wall, hSpself, hSpcross, hDpself, hDpcross);

  % 4b. CALDERON right-preconditioner: RANK-COMPLETED on-surface S value apply,
  %     P = S + nh1 nh1' + nh2 nh2' (bare Stokes S has the pressure null space S[n_k]=0, so
  %     sigma = S tau could not carry the normal components; the completion restores them).
  %     D'S = K'^2 - 1/4 I -> (-1/2+S'+D')P is second-kind + low-rank -> iterations flat.
  hSvself = axpso_close_setup_mex(2, 1, 1, 0, mu, 1, K2, pv, npv, pmv, iside, iclosed, gate, ...
              geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
              rball, K2*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);            % SLP value self
  axpso_close2corr_set_mex(hSvself, K2, geomoff, nsx, sxf, snxf, swsf, targoff, K2*Nnod, Xall, Nall, [C{:}]);
  hSvcross = axpso_close_setup_mex(2, 1, 1, 0, mu, 2, K2, pv, npv, pmv, iside, iclosed, gate, ...
               geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
               rball, K2*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);           % SLP value cross
  axpso_close2corr_set_mex(hSvcross, K2, geomoff, nsx, sxf, snxf, swsf, targoff, K2*Nnod, Xall, Nall, [C{:}]);
  nh1=zeros(3*Mtot,1); nh1(blk{1})=NX{1}(:); nh1=nh1/norm(nh1);
  nh2=zeros(3*Mtot,1); nh2(blk{2})=NX{2}(:); nh2=nh2/norm(nh2);
  applyP=@(x) applymvS(x, fmmeps, Xall, wall, hSvself, hSvcross) + nh1*(nh1'*x) + nh2*(nh2'*x);

  % 5. solve (gmres on the SECOND-KIND composite (-1/2+S'+D')P, then sigma = P tau)
  [tau,flag,relres,iter]=gmres(@(t) applyA(applyP(t)), g_tr, [], 1e-11, 300);
  sigma=applyP(tau);
  fprintf('  gmres: flag=%d  iters=%d  relres=%.3e\n', flag, iter(2), relres);

  % 6. 3D exterior target grid (outside both particles)
  cmid=(C{1}+C{2})/2; gv=linspace(-3,3,22); [GX,GY,GZ]=ndgrid(gv,gv,gv); P=[GX(:)';GY(:)';GZ(:)']+cmid;
  ext=true(1,size(P,2)); dmin=inf(1,size(P,2));
  for k=1:2, loc=R{k}.'*(P-C{k}); sd=(hypot(loc(1,:),loc(2,:))/ac(k,1)).^2+(loc(3,:)/ac(k,2)).^2; ext=ext&(sd>1.05); dmin=min(dmin,sd-1); end
  Pe=P(:,ext); Me=size(Pe,2); dmin=dmin(ext);

  % 7. field eval u = (S+D)[sigma]: FMM naive (S+D) at Pe + {S,D} eval corrections (targoff=ones ->
  %    every source sees ALL Me grid targets), same full-minus-naive recipe applied source 3*K*Nnod -> 3*Me.
  u = Sto3dSLPfmm_il(struct('x',Pe), struct('x',Xall,'w',wall), sigma, fmmeps) ...
    + Sto3dDLPfmm_il(struct('x',Pe), struct('x',Xall,'w',wall,'nx',Nall), sigma, fmmeps);
  hSeval = axpso_close_setup_mex(2, 1, 1, 0, mu, 2, K2, pv, npv, pmv, iside, iclosed, gate, ...
             geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
             rball, Me, ones(K2+1,1), Pe, Pe, cat(3,R{:}), [C{:}], 0);                % SLP eval
  axpso_close2corr_set_mex(hSeval, K2, geomoff, nsx, sxf, snxf, swsf, ones(K2+1,1), Me, Pe, Pe, [C{:}]);
  hDeval = axpso_close_setup_mex(2, 3, 1, 0, mu, 2, K2, pv, npv, pmv, iside, iclosed, gate, ...
             geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
             rball, Me, ones(K2+1,1), Pe, Pe, cat(3,R{:}), [C{:}], 0);                % DLP eval
  axpso_close2corr_set_mex(hDeval, K2, geomoff, nsx, sxf, snxf, swsf, ones(K2+1,1), Me, Pe, Pe, [C{:}]);
  u=axpso_corr_apply_mex(hSeval, 3*K2*Nnod, sigma, 3*Me, u);
  u=axpso_corr_apply_mex(hDeval, 3*K2*Nnod, sigma, 3*Me, u);
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
  xlabel('x'); ylabel('y'); zlabel('z'); title('2-particle Stokes DLPn (S''+D'', physmat sparse-correction): 3D velocity error (last refinement)');
  figure(2),ylim([-1 3])
end

figure(2),clf; hold on;
scatter3(Pe(1,:),Pe(2,:),Pe(3,:),14,log10(max(err,1e-17)),'filled');
plot3(X{1}(1,:),X{1}(2,:),X{1}(3,:),'.','Color',[.6 .6 .6]);
plot3(X{2}(1,:),X{2}(2,:),X{2}(3,:),'.','Color',[.6 .6 .6]);
plot3(Y(1,:),Y(2,:),Y(3,:),'p','MarkerSize',16,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
cb=colorbar; cb.Label.String='log_{10}|u_h-u_{exact}|/max|u_{exact}|';
clim([-16 -4]); axis equal; view(35,18); grid on; colormap('jet');
xlabel('x'); ylabel('y'); zlabel('z'); title('2-particle Stokes DLPn (S''+D'', physmat sparse-correction): 3D velocity error (last refinement)');
figure(2),ylim([-1 3])

% --- convergence plot: max field error vs panel count ---
figure(1),clf; semilogy(np_vals,errmax,'o-k'); grid on; xlabel('n_p'); ylabel('max field err');
title('multi-particle Stokes DLPn (S''+D'', physmat sparse-correction): h-refinement, p=16, M=2 n_p');

% cmp = getPyPlot_cMap('rainbow', [], [], '"/Users/hzhu/.pyenv/versions/3.11.13/bin/python"');
% colormap(cmp)

% exportgraphics(figure(1),'axissymsstok_stok_dlpn_multi_physmat_convergence.png','Resolution',200)
% exportgraphics(figure(2),'axissymsstok_stok_dlpn_multi_physmat_error.png','Resolution',200)

function y = applymvSD(x, fmmeps, Xall, Nall, wall, hSpself, hSpcross, hDpself, hDpcross)
% matrix-free (-1/2 I + S' + D') matvec: ONE stfmm3d call (stoklet + strslet/strsvec) -> traction of
% (S+D) at the sources t = (grad u + grad u^T) n - p n (self i~=j excluded; the -1/2 jump rides in
% the S' self correction) + four handle corr applies (ACCUMULATE, R sandwich in-module).
sig3=reshape(x,3,[]).*wall;                                   % weighted density (Stokeslet + stresslet strength)
srcinfo=[]; srcinfo.sources=Xall; srcinfo.stoklet=sig3; srcinfo.strslet=sig3; srcinfo.strsvec=Nall;
U=stfmm3d(fmmeps,srcinfo,3);                                  % pot + pre + grad at sources
Gs=U.grad+permute(U.grad,[2 1 3]);                            % only the symmetric part enters
y=squeeze(sum(Gs.*reshape(Nall,1,3,[]),2)) - U.pre.*Nall;
y=y(:);
n=numel(x);
y=axpso_corr_apply_mex(hSpself,  n, x, n, y);
y=axpso_corr_apply_mex(hSpcross, n, x, n, y);
y=axpso_corr_apply_mex(hDpself,  n, x, n, y);
y=axpso_corr_apply_mex(hDpcross, n, x, n, y);
end

function y = applymvS(x, fmmeps, Xall, wall, hSvself, hSvcross)
% on-surface S VALUE apply (Calderon preconditioner core): Sto3dSLPfmm_il naive velocity at the
% sources (targets==sources -> self dropped) + {self,cross} corr applies.
y = Sto3dSLPfmm_il(struct('x',Xall), struct('x',Xall,'w',wall), x, fmmeps);
n=numel(x);
y=axpso_corr_apply_mex(hSvself,  n, x, n, y);
y=axpso_corr_apply_mex(hSvcross, n, x, n, y);
end
