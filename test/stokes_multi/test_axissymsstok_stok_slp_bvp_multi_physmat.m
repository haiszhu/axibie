clearvars; format short e;
addpath('../../../axibie/utils');
addpath('../../utils');
addpath('../../matlab');
addpath('../../external/fmm3d/matlab');                       % stfmm3d (via Sto3dSLPfmm_il)

% MULTI-particle Stokes SLP exterior-Dirichlet BVP -- PHYSICAL-SPACE SPARSE-CORRECTION version of
% test_axissymsstok_stok_slp_bvp_multi.m, in the framework style of the SLPn/DLP physmat twins.
% Single-layer representation u = S[sigma]; match the surface velocity S[sigma] = u_ex.  Fully
% MATRIX-FREE mex; each near correction is FULL close-eval minus naive, in Fortran:
%   SETUP (per pass): axpso_close_setup_mex returns the FULL close-eval matrix (iform=0) via the
%     WHOLE-PARTICLE split-open closestokslp build, then axpso_close2corr_set_mex subtracts the naive
%     Sto3d Stokeslet IN FORTRAN (compact layout, self node auto-zero) -> the handle holds the
%     correction.  Two SYSTEM handles: {self, cross}.
%   SYSTEM (velocity): Sto3dSLPfmm_il naive S at the sources (targets==sources -> self dropped) +
%     {self,cross} corr applies (ACCUMULATE, lab<->local R sandwich INSIDE the module).
%   SYSTEM (explicit): naive dense Sto3dSLPmat_il (self diag zeroed) + the handle corrections read
%     back via axpso_corr_get_mex and scattered exactly as corr_apply would (self: circulant +
%     R(phi_a) conjugation; cross: global target rows via corr_get's canon->global map).
%   SOLVE: DIRECT lsqminnorm (S is first-kind AND rank-deficient -- per-particle pressure null
%     space -- so min-norm, == the dense reference; no gmres, no deflation).
%   EVAL (velocity): Sto3dSLPfmm_il naive + SLP eval correction, same recipe.
% VECTOR (nc=3), node-interleaved lab frame; mu=1.

p=16; np_vals=2:8; errmax=nan(1,numel(np_vals)); mu=1.0; gate=2.0; fmmeps=1e-12;   % EXPLICIT dense (3*2*Nnod)^2: np=6 ~ 15 GB
rot=@(u,th) cos(th)*eye(3)+sin(th)*[0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0]+(1-cos(th))*(u*u');
Zf={@(t) 0.5*sin(t)-1i*1.0*cos(t), @(t) 1.0*sin(t)-1i*0.5*cos(t)}; ac=[0.5 1.0; 1.0 0.5];
stk=@(P,y,Fv)(1/(8*pi*mu))*( Fv./vecnorm(P-y) + (Fv.'*(P-y)).*(P-y)./vecnorm(P-y).^3 );    % interior Stokeslet velocity
yloc={[0.10 0.05; 0 0.05; 0.20 -0.30],[0.20 -0.10; 0 0.15; 0.10 -0.10]};
Floc={[1.0 -0.5; -0.7 0.6; 0.5 -0.4],[0.6 -0.4; 0.3 -0.5; -0.7 0.5]};

for kk=1:numel(np_vals)
  np=np_vals(kk); M=2*np; nphi=2*M+1; ph=2*pi*(0:nphi-1)/nphi; iside=1; iclosed=0;
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
  rnear=gate*sqrt(max(arrayfun(@(j) sum(geo{1}.sws((j-1)*p+(1:p))), 1:np)));

  % 3. SYSTEM corrections: FULL close-eval (close_setup ikernel=2 ilayer=1, whole-particle build) minus
  %    the naive Sto3d Stokeslet (close2corr_set, in Fortran) -- self and cross.
  hSself = axpso_close_setup_mex(2, 1, 1, 0, mu, 1, K2, pv, npv, pmv, iside, iclosed, gate, ...
             geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
             1.0+rnear, K2*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);         % SLP self
  axpso_close2corr_set_mex(hSself, K2, geomoff, nsx, sxf, snxf, swsf, targoff, K2*Nnod, Xall, Nall, [C{:}]);
  hScross = axpso_close_setup_mex(2, 1, 1, 0, mu, 2, K2, pv, npv, pmv, iside, iclosed, gate, ...
              geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
              1.0+rnear, K2*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);        % SLP cross
  axpso_close2corr_set_mex(hScross, K2, geomoff, nsx, sxf, snxf, swsf, targoff, K2*Nnod, Xall, Nall, [C{:}]);

  % 4. EXPLICIT dense operator: naive Sto3dSLPmat_il (self 3x3 node-diag zeroed) + the handle
  %    corrections scattered IN FORTRAN by axpso_corr2dense_get_mex (the explicit-matrix twin of
  %    corr_apply: self = circulant + R(phi_a) conjugation + R_k sandwich; cross = global canon rows).
  A = real(Sto3dSLPmat_il(struct('x',Xall), struct('x',Xall,'w',wall)));   % naive full S velocity
  A = zeroselfblk(A, Mtot);
  A = axpso_corr2dense_get_mex(hSself,  3*Mtot, 3*Mtot, A);
  A = axpso_corr2dense_get_mex(hScross, 3*Mtot, 3*Mtot, A);

  % 5. solve (DIRECT: S is rank-deficient -- per-particle pressure null space -> min-norm)
  sigma = lsqminnorm(A, g_di);
  fprintf('  lsqminnorm: |A sigma - g| = %.3e\n', norm(A*sigma-g_di)/norm(g_di));

  % 6. 3D exterior target grid (outside both particles)
  cmid=(C{1}+C{2})/2; gv=linspace(-3,3,22); [GX,GY,GZ]=ndgrid(gv,gv,gv); P=[GX(:)';GY(:)';GZ(:)']+cmid;
  ext=true(1,size(P,2)); dmin=inf(1,size(P,2));
  for k=1:2, loc=R{k}.'*(P-C{k}); sd=(hypot(loc(1,:),loc(2,:))/ac(k,1)).^2+(loc(3,:)/ac(k,2)).^2; ext=ext&(sd>1.05); dmin=min(dmin,sd-1); end
  Pe=P(:,ext); Me=size(Pe,2); dmin=dmin(ext);

  % 7. field eval u = S[sigma]: Sto3dSLPfmm_il naive velocity + SLP eval correction (targoff=ones ->
  %    every source sees ALL Me grid targets), applied source 3*K*Nnod -> 3*Me (R sandwich in-module).
  u = Sto3dSLPfmm_il(struct('x',Pe), struct('x',Xall,'w',wall), sigma, fmmeps);
  hSeval = axpso_close_setup_mex(2, 1, 1, 0, mu, 2, K2, pv, npv, pmv, iside, iclosed, gate, ...
             geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
             1.0+rnear, Me, ones(K2+1,1), Pe, Pe, cat(3,R{:}), [C{:}], 0);
  axpso_close2corr_set_mex(hSeval, K2, geomoff, nsx, sxf, snxf, swsf, ones(K2+1,1), Me, Pe, Pe, [C{:}]);
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
  xlabel('x'); ylabel('y'); zlabel('z'); title('2-particle Stokes SLP (physmat sparse-correction): 3D velocity error (last refinement)');
  figure(2),ylim([-1 3])
end

figure(2),clf; hold on;
scatter3(Pe(1,:),Pe(2,:),Pe(3,:),14,log10(max(err,1e-17)),'filled');
plot3(X{1}(1,:),X{1}(2,:),X{1}(3,:),'.','Color',[.6 .6 .6]);
plot3(X{2}(1,:),X{2}(2,:),X{2}(3,:),'.','Color',[.6 .6 .6]);
plot3(Y(1,:),Y(2,:),Y(3,:),'p','MarkerSize',16,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
cb=colorbar; cb.Label.String='log_{10}|u_h-u_{exact}|/max|u_{exact}|';
clim([-16 -4]); axis equal; view(35,18); grid on; colormap('jet');
xlabel('x'); ylabel('y'); zlabel('z'); title('2-particle Stokes SLP (physmat sparse-correction): 3D velocity error (last refinement)');
figure(2),ylim([-1 3])

% --- convergence plot: max field error vs panel count ---
figure(1),clf; semilogy(np_vals,errmax,'o-k'); grid on; xlabel('n_p'); ylabel('max field err');
title('multi-particle Stokes SLP (physmat sparse-correction): h-refinement, p=16, M=2 n_p');

% cmp = getPyPlot_cMap('rainbow', [], [], '"/Users/hzhu/.pyenv/versions/3.11.13/bin/python"');
% colormap(cmp)

% exportgraphics(figure(1),'axissymsstok_stok_slp_multi_physmat_convergence.png','Resolution',200)
% exportgraphics(figure(2),'axissymsstok_stok_slp_multi_physmat_error.png','Resolution',200)

function A = zeroselfblk(A, m)
% zero the m node-diagonal 3x3 blocks of a 3m x 3m node-interleaved matrix (singular self-interaction)
N3 = 3*m; ri = (1:3).' + 3*(0:m-1);
lin = (reshape(ri,3,1,m)-1)*N3 + reshape(ri,1,3,m);
A(lin(:)) = 0;
end
