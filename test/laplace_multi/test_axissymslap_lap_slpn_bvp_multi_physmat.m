clearvars; format short e;
addpath('../../../axibie/utils');
addpath('../../utils');
addpath('../../matlab');
addpath('../../external/fmm3d/matlab');                       % lfmm3d / Lap3dSLPfmm

% MULTI-particle Laplace SLPn (S') exterior-Neumann BVP -- PHYSICAL-SPACE SPARSE-CORRECTION version
% of test_axissymslap_lap_slpn_bvp_multi.m, in the framework style of the Stokes twin
% test_axissymsstok_stok_slpn_bvp_multi_physmat.m.  Same setup (2 random-pose ellipsoids, interior
% point charges), but the operator is fully MATRIX-FREE mex and the near correction is built via the
% FULL close-eval operator minus the naive kernel:
%   SETUP (per pass): axpso_close_setup_mex returns the FULL Laplace close-eval matrix (iform=0);
%     axpso_close2corr_set_mex subtracts the naive Lap3d kernel IN FORTRAN (compact layout, self node
%     auto-zero) -> the handle now holds the correction == what corr_setup would store.
%   SYSTEM (flux): lfmm3d naive S' = grad(pot).n at the sources (self i~=j excluded by the FMM; the
%     -1/2 exterior jump rides in the SELF correction) + hself/hcross corr apply (each ACCUMULATES).
%   SOLVE: gmres, NO deflation -- (-1/2 + S') is second-kind FULL-RANK for Laplace exterior-Neumann.
%   EVAL (potential): Lap3dSLPfmm naive + heval (ilayer=1 SLP value) correction, same recipe.
% SCALAR (nc=1): one dof/node, no lab<->local R sandwich.

p=16; np_vals=2:12; errmax=nan(1,numel(np_vals)); gate=2.0; fmmeps=1e-12;
rot=@(u,th) cos(th)*eye(3)+sin(th)*[0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0]+(1-cos(th))*(u*u');
Zf={@(t) 0.5*sin(t)-1i*1.0*cos(t), @(t) 1.0*sin(t)-1i*0.5*cos(t)}; ac=[0.5 1.0; 1.0 0.5];
yloc={[0.10 0.05; 0 0.05; 0.20 -0.30],[0.20 -0.10; 0 0.15; 0.10 -0.10]}; qq=[1.0 -0.7 0.6 -0.5];

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
  N=numel(geo{1}.sx); Nnod=N*nphi; Mtot=2*Nnod;
  Xall=[X{1} X{2}]; Nall=[NX{1} NX{2}]; wall=[sq{1}.w sq{2}.w];
  nodeblk={1:Nnod, Nnod+1:2*Nnod};

  % interior point charges -> exact exterior potential uex (eval) AND exact surface flux dn(uex) (Neumann data)
  Y=zeros(3,0); for k=1:2, Y=[Y, R{k}*yloc{k}+C{k}]; end
  uex=@(P) lapsum(P,Y,qq);
  g_tr=zeros(1,Mtot);                                          % dn(uex) = -sum_j q_j (n.(X-Y_j))/(4pi r^3)
  for j=1:numel(qq), dXY=Xall-Y(:,j); rr=vecnorm(dXY); g_tr=g_tr - qq(j)*(sum(Nall.*dXY,1))./(4*pi*rr.^3); end
  g_tr=g_tr(:);

  % 2. flat per-particle meridian geometry consumed by the handle-based setups (self + cross)
  K2=2;
  sxs=zeros(N,K2); snxs=sxs; swss=zeros(N,K2); swxps=zeros(N,K2); tpans=zeros(np+1,K2);
  for k=1:K2, gs=geo{k}; sxs(:,k)=gs.sx; snxs(:,k)=gs.snx; swss(:,k)=gs.sws; swxps(:,k)=gs.swxp; tpans(:,k)=gs.tpan; end
  pv=p*ones(K2,1); npv=np*ones(K2,1); pmv=M*ones(K2,1);
  geomoff=1+(0:K2)'*N; tpanoff=1+(0:K2)'*(np+1); targoff=1+(0:K2)'*Nnod;
  nsx=K2*N; ntpan=K2*(np+1);
  sxf=sxs(:); snxf=snxs(:); swsf=swss(:); swxpf=swxps(:); tpanf=tpans(:);
  rnear=gate*sqrt(max(arrayfun(@(j) sum(geo{1}.sws((j-1)*p+(1:p))), 1:np)));

  % 3. SELF corrections: FULL close-eval (close_setup ikernel=1 laplace, ilayer=2 SLPn) minus the naive
  %    Lap3d S' flux (close2corr_set, in Fortran; the self node auto-zeros).  hself holds the correction.
  hself = axpso_close_setup_mex(1, 2, 1, 0, 0.0, 1, K2, pv, npv, pmv, iside, iclosed, gate, ...
            geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
            1.0+rnear, K2*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);
  axpso_close2corr_set_mex(hself, K2, geomoff, nsx, sxf, snxf, swsf, targoff, K2*Nnod, Xall, Nall, [C{:}]);
  % peek: read a stored block back per particle (Sk = the near CORRECTION matrix, full - naive)
  for k=1:K2
    gs=geo{k}; t3d0=[real(gs.sx).'; zeros(1,N); imag(gs.sx).'];
    tcxi0=zeros(np+1,1); nt=0;
    [~,nt]=axps_closesize_mex(N, gs.sx, t3d0, p, np, gs.sx, gs.sws, gate, tcxi0, nt);
    [bg, ig, tg] = axpso_corr_get_mex(hself, k, nt, 1, p, np, M, zeros(nt*nphi*p,1), zeros(nt,1), zeros(np+1,1));
    Sk = reshape(bg, nt, nphi*p);
    fprintf('  self block k=%d:  %d x %d  (ntcx=%d)\n', k, size(Sk,1), size(Sk,2), nt);
  end

  % 4. CROSS corrections: FULL close-eval (iinter=2, targets = the other particle's nodes, own block
  %    dropped via targoff, rotated into each source frame in-module) minus the naive S' flux.
  hcross = axpso_close_setup_mex(1, 2, 1, 0, 0.0, 2, K2, pv, npv, pmv, iside, iclosed, gate, ...
             geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
             1.0+rnear, K2*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);
  axpso_close2corr_set_mex(hcross, K2, geomoff, nsx, sxf, snxf, swsf, targoff, K2*Nnod, Xall, Nall, [C{:}]);

  % 5. matrix-free operator: lfmm3d naive S' + hself/hcross corr apply (scalar, no deflation)
  applyA=@(x) applymv(x, fmmeps, Xall, Nall, wall, hself, hcross);

  % 6. solve (gmres, second-kind full-rank -> no deflation)
  [sigma,flag,relres,iter]=gmres(applyA, g_tr, [], 1e-14, 300);
  fprintf('  gmres: flag=%d  iters=%d  relres=%.3e\n', flag, iter(2), relres);

  % 7. 3D exterior target grid (outside both particles)
  cmid=(C{1}+C{2})/2; gv=linspace(-3,3,22); [GX,GY,GZ]=ndgrid(gv,gv,gv); P=[GX(:)';GY(:)';GZ(:)']+cmid;
  ext=true(1,size(P,2)); dmin=inf(1,size(P,2));
  for k=1:2, loc=R{k}.'*(P-C{k}); sd=(hypot(loc(1,:),loc(2,:))/ac(k,1)).^2+(loc(3,:)/ac(k,2)).^2; ext=ext&(sd>1.05); dmin=min(dmin,sd-1); end
  Pe=P(:,ext); Me=size(Pe,2); dmin=dmin(ext);

  % 8. field eval u = S[sigma]: Lap3dSLPfmm naive potential + heval (ilayer=1 SLP value) correction,
  %    FULL close-eval minus naive SLP potential (close2corr_set), applied source K*Nnod -> target Me.
  u=Lap3dSLPfmm(struct('x',Pe), struct('x',Xall,'w',wall), sigma, fmmeps);
  heval = axpso_close_setup_mex(1, 1, 1, 0, 0.0, 2, K2, pv, npv, pmv, iside, iclosed, gate, ...
            geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
            1.0+rnear, Me, ones(K2+1,1), Pe, Pe, cat(3,R{:}), [C{:}], 0);
  axpso_close2corr_set_mex(heval, K2, geomoff, nsx, sxf, snxf, swsf, ones(K2+1,1), Me, Pe, Pe, [C{:}]);
  u=axpso_corr_apply_mex(heval, K2*Nnod, sigma, Me, u);
  Uex=uex(Pe).'; err=abs(u-Uex)/max(abs(Uex)); inb=dmin.'<0.1; errmax(kk)=max(err);
  fprintf('np=%2d  N=%5d/part  M=%2d  M3=%5d :  near(d<0.1) %.2e   far(d>=0.1) %.2e\n', ...
          np, Nnod, M, Me, max([err(inb);0]), max(err(~inb)));

  figure(2),clf; hold on;
  scatter3(Pe(1,:),Pe(2,:),Pe(3,:),14,log10(max(err,1e-17)),'filled');
  plot3(X{1}(1,:),X{1}(2,:),X{1}(3,:),'.','Color',[.6 .6 .6]);
  plot3(X{2}(1,:),X{2}(2,:),X{2}(3,:),'.','Color',[.6 .6 .6]);
  plot3(Y(1,:),Y(2,:),Y(3,:),'p','MarkerSize',16,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
  cb=colorbar; cb.Label.String='log_{10}|u_h-u_{exact}|/max|u_{exact}|';
  axis equal; view(35,18); grid on; colormap('jet');
  xlabel('x'); ylabel('y'); zlabel('z'); title('2-particle Laplace SLPn (S'', physmat sparse-correction): 3D field error (last refinement)');
  figure(2),ylim([-1 3])

end

figure(2),clf; hold on;
scatter3(Pe(1,:),Pe(2,:),Pe(3,:),14,log10(max(err,1e-17)),'filled');
plot3(X{1}(1,:),X{1}(2,:),X{1}(3,:),'.','Color',[.6 .6 .6]);
plot3(X{2}(1,:),X{2}(2,:),X{2}(3,:),'.','Color',[.6 .6 .6]);
plot3(Y(1,:),Y(2,:),Y(3,:),'p','MarkerSize',16,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
cb=colorbar; cb.Label.String='log_{10}|u_h-u_{exact}|/max|u_{exact}|';
clim([-16 -4]); axis equal; view(35,18); grid on; colormap('jet');
xlabel('x'); ylabel('y'); zlabel('z'); title('2-particle Laplace SLPn (S'', physmat sparse-correction): 3D field error (last refinement)');
figure(2),ylim([-1 3])

% --- convergence plot: max field error vs panel count ---
figure(1),clf; semilogy(np_vals,errmax,'o-k'); grid on; xlabel('n_p'); ylabel('max field err');
title('multi-particle Laplace SLPn (S'', physmat sparse-correction): h-refinement, p=16, M=2 n_p');

% cmp = getPyPlot_cMap('rainbow', [], [], '"/Users/hzhu/.pyenv/versions/3.11.13/bin/python"');
% colormap(cmp)

% exportgraphics(figure(1),'axissymslap_lap_slpn_multi_physmat_convergence.png','Resolution',200)
% exportgraphics(figure(2),'axissymslap_lap_slpn_multi_physmat_error.png','Resolution',200)

function y = applymv(x, fmmeps, Xall, Nall, wall, hself, hcross)
% matrix-free (-1/2 I + S') matvec: global lfmm3d naive normal-derivative (self i~=j excluded; the
% -1/2 jump rides in the self correction) + hself/hcross corr apply (each ACCUMULATES).  SCALAR.
srcinfo.sources=Xall; srcinfo.charges=(x(:).').*wall;         % scalar charges, weights baked in
U=lfmm3d(fmmeps,srcinfo,2);                                   % pot + grad at sources (self i~=j excluded)
y=(sum(U.grad.*Nall,1)).';                                    % S' = grad . node normal, scalar per node
n=numel(x);
y=axpso_corr_apply_mex(hself,  n, x, n, y);                  % self correction
y=axpso_corr_apply_mex(hcross, n, x, n, y);                  % cross correction
end

function u=lapsum(P,Y,Q)
u=zeros(1,size(P,2));
for j=1:size(Y,2), u=u+Q(j)./(4*pi*vecnorm(P-Y(:,j))); end
end
