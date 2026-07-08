clearvars; format short e;
addpath('../../utils');
addpath('../../matlab');
addpath('../../external/fmm3d/matlab');                       % stfmm3d / Sto3dSLPfmm

p=16; np_vals=2:8; errmax=nan(1,numel(np_vals)); mu=1.0; gate=2.0;
rot=@(u,th) cos(th)*eye(3)+sin(th)*[0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0]+(1-cos(th))*(u*u');
Zf={@(t) 0.5*sin(t)-1i*1.0*cos(t), @(t) 1.0*sin(t)-1i*0.5*cos(t)}; ac=[0.5 1.0; 1.0 0.5];
stk=@(P,y,Fv)(1/(8*pi))*( Fv./vecnorm(P-y) + (Fv.'*(P-y)).*(P-y)./vecnorm(P-y).^3 );    % interior Stokeslet velocity
yloc={[0.10 0.05; 0 0.05; 0.20 -0.30],[0.20 -0.10; 0 0.15; 0.10 -0.10]};
Floc={[1.0 -0.5; -0.7 0.6; 0.5 -0.4],[0.6 -0.4; 0.3 -0.5; -0.7 0.5]};

for kk=1:numel(np_vals)
  np=np_vals(kk); M=2*np; nphi=2*M+1; ph=2*pi*(0:nphi-1)/nphi; iside=1; iclosed=0;
  rng(1);

  % 1. particles: random pose; LOCAL-frame meridian geometry + 3D ring nodes/normals (lab + local)
  C=cell(2,1); R=cell(2,1); X=cell(2,1); NX=cell(2,1); sq=cell(2,1); geo=cell(2,1); snx3=cell(2,1);
  for k=1:2
    s=[]; s.p=p; s.Z=Zf{k}; s.tpan=linspace(0,pi,np+1)'; s=quadr(s,[],'p','G');
    d=randn(3,1); C{k}=(k-1)*1.6*d/norm(d); uu=randn(3,1); uu=uu/norm(uu); R{k}=rot(uu,2*pi*rand);
    rho=real(s.x); z=imag(s.x); nr=real(s.nx); nz=imag(s.nx);
    X{k}=R{k}*[(kron(cos(ph(:)),rho)).';(kron(sin(ph(:)),rho)).';(kron(ones(nphi,1),z)).']+C{k};
    NX{k}=R{k}*[(kron(cos(ph(:)),nr)).';(kron(sin(ph(:)),nr)).';(kron(ones(nphi,1),nz)).'];   % 3D node normals (lab)
    sq{k}=axisym_to_3d_quadrature(real(s.x),imag(s.x),s.ws,nphi);                             % LOCAL-frame 3D nodes/weights
    snx3{k}=[reshape(real(s.nx(:))*cos(ph),1,[]); reshape(real(s.nx(:))*sin(ph),1,[]); repmat(imag(s.nx(:)).',1,nphi)];
    g=[]; g.sx=s.x(:); g.snx=s.nx(:); g.sws=s.ws(:); g.swxp=s.wxp(:); g.tpan=s.tpan(:); geo{k}=g;
  end
  N=numel(geo{1}.sx); Nnod=N*nphi; Ndof=3*Nnod; Mtot=2*Nnod;
  Xall=[X{1} X{2}]; Nall=[NX{1} NX{2}]; wall=[sq{1}.w sq{2}.w];
  nodeblk={1:Nnod, Nnod+1:2*Nnod}; blk={1:Ndof, Ndof+1:2*Ndof};

  % interior Stokeslets -> exact exterior velocity uex (for eval) AND exact surface traction (Neumann data)
  Y=zeros(3,0); Fg=zeros(3,0); for k=1:2, Y=[Y, R{k}*yloc{k}+C{k}]; Fg=[Fg, R{k}*Floc{k}]; end
  uex=@(P) stk(P,Y(:,1),Fg(:,1))+stk(P,Y(:,2),Fg(:,2))+stk(P,Y(:,3),Fg(:,3))+stk(P,Y(:,4),Fg(:,4));
  [~,SpY]=Sto3dSLPmat_il(struct('x',Xall,'nx',Nall), struct('x',Y,'w',ones(1,size(Y,2))));
  g = SpY*Fg(:);                                              % exact surface traction, node-interleaved (lab)

  % 2. SELF corrections, formed EXPLICITLY in this driver from the FULL close-eval operator:
  %    axpso_close_setup_mex returns the FULL close-eval matrix entries (iform=0), NOT the
  %    correction.  close2corr then reads each block (corr_get), subtracts the naive Sto3d SLPn
  %    traction it computes here (the self 3x3 node-diagonal auto-zeros at the coincident node),
  %    and writes the difference back (corr_set) so the handle now holds the correction = full -
  %    naive == what corr_setup would have stored.  Thereafter corr_apply scatters it identically.
  K2=2;
  sxs=zeros(N,K2); snxs=sxs; swss=zeros(N,K2); swxps=zeros(N,K2); tpans=zeros(np+1,K2);
  for k=1:K2, gs=geo{k}; sxs(:,k)=gs.sx; snxs(:,k)=gs.snx; swss(:,k)=gs.sws; swxps(:,k)=gs.swxp; tpans(:,k)=gs.tpan; end
  pv=p*ones(K2,1); npv=np*ones(K2,1); pmv=M*ones(K2,1);
  geomoff=1+(0:K2)'*N; tpanoff=1+(0:K2)'*(np+1); targoff=1+(0:K2)'*Nnod;
  nsx=K2*N; ntpan=K2*(np+1);
  sxf=sxs(:); snxf=snxs(:); swsf=swss(:); swxpf=swxps(:); tpanf=tpans(:);
  rnear=gate*sqrt(max(arrayfun(@(j) sum(geo{1}.sws((j-1)*p+(1:p))), 1:np)));
  [hself,~]=axpso_close_setup_mex(2, 2, mu, 1, K2, pv, npv, pmv, iside, iclosed, gate, ...
            geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
            1.0+rnear, K2*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);
  axpso_close2corr_set_mex(hself, K2, geomoff, nsx, sxf, snxf, swsf, targoff, K2*Nnod, Xall, Nall, [C{:}]);   % full - naive -> correction (Fortran)
  % peek: read a stored block back per particle (Sk = the near CORRECTION matrix, full - naive).
  for k=1:K2
    gs=geo{k}; t3d0=[real(gs.sx).'; zeros(1,N); imag(gs.sx).'];
    tcxi0=zeros(np+1,1); nt=0;
    [~,nt]=axps_closesize_mex(N, gs.sx, t3d0, p, np, gs.sx, gs.sws, gate, tcxi0, nt);
    [bg, ig, tg] = axpso_corr_get_mex(hself, k, nt, 3, p, np, M, ...
          zeros(3*3*nt*nphi*p,1), zeros(nt,1), zeros(np+1,1));
    Sk = reshape(bg, 3*nt, 3*nphi*p);
    fprintf('  self block k=%d:  %d x %d  (ntcx=%d)\n', k, size(Sk,1), size(Sk,2), nt);
  end

  % 3. CROSS corrections, same recipe: FULL close-eval (iinter=2, targets = the other particle's
  %    nodes, own block dropped via targoff, rotated into each source frame inside the module) minus
  %    the naive SLPn traction, formed here (corr_get returns each near row's GLOBAL target index).
  [hcross,~]=axpso_close_setup_mex(2, 2, mu, 2, K2, pv, npv, pmv, iside, iclosed, gate, ...
             geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
             1.0+rnear, K2*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);
  axpso_close2corr_set_mex(hcross, K2, geomoff, nsx, sxf, snxf, swsf, targoff, K2*Nnod, Xall, Nall, [C{:}]);  % full - naive (Fortran)

  % 4. matrix-free operator: stfmm3d naive S' + hself/hcross corr apply (R sandwich INSIDE the
  %    module) + per-particle pressure-gauge deflation (S[n_k]=0, cf the single-particle physmat2)
  nh1=zeros(3*Mtot,1); nh1(blk{1})=NX{1}(:); nh1=nh1/norm(nh1);
  nh2=zeros(3*Mtot,1); nh2(blk{2})=NX{2}(:); nh2=nh2/norm(nh2);
  applyA=@(x) applymv(x, Xall, Nall, wall, hself, hcross) + nh1*(nh1'*x) + nh2*(nh2'*x);

  % 5. solve (unrestarted GMRES on the deflated matrix-free operator)
  [sigma,flag,relres,iter]=gmres(applyA, g, [], 1e-11, 300);
  fprintf('  gmres: flag=%d  iters=%d  relres=%.3e  sig.n=[%.2e %.2e]\n', ...
          flag, iter(2), relres, nh1'*sigma, nh2'*sigma);

  % 6. 3D exterior target grid (outside both particles) -- same as the dense multi test
  cmid=(C{1}+C{2})/2; gv=linspace(-3,3,22); [GX,GY,GZ]=ndgrid(gv,gv,gv); P=[GX(:)';GY(:)';GZ(:)']+cmid;
  ext=true(1,size(P,2)); dmin=inf(1,size(P,2));
  for k=1:2, loc=R{k}.'*(P-C{k}); sd=(hypot(loc(1,:),loc(2,:))/ac(k,1)).^2+(loc(3,:)/ac(k,2)).^2; ext=ext&(sd>1.05); dmin=min(dmin,sd-1); end
  Pe=P(:,ext); Me=size(Pe,2); dmin=dmin(ext);

  % 7. field eval u = S[sigma]: Sto3dSLPfmm naive velocity seed + eval correction, formed here by
  %    the SAME full-minus-naive recipe: heval holds the FULL SLP-velocity close-eval (close_setup
  %    ilayer=1, targoff=ones -> every source sees ALL Me targets, rotated into its frame inside the
  %    module); close2corr subtracts the naive Sto3d SLP velocity (ilayer=1, no target normals) and
  %    writes back the correction, then corr_apply adds it to the FMM seed (R sandwich in-module).
  sigblk = reshape(reshape(sigma,3,[]).',[],1);               % interlocked -> component-block
  ublk   = Sto3dSLPfmm(struct('x',Pe), struct('x',Xall,'w',wall), sigblk, 1e-14);
  ucm    = reshape(reshape(ublk,[],3).',[],1);                % component-block -> interlocked (3*Me)
  heval = axpso_close_setup_mex(2, 1, mu, 2, K2, pv, npv, pmv, iside, iclosed, gate, ...
            geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
            1.0+rnear, Me, ones(K2+1,1), Pe, Pe, cat(3,R{:}), [C{:}], 0);
  axpso_close2corr_set_mex(heval, K2, geomoff, nsx, sxf, snxf, swsf, ones(K2+1,1), Me, Pe, Pe, [C{:}]);   % full - naive (SLP, Fortran)
  ucm = axpso_corr_apply_mex(heval, 3*K2*Nnod, sigma, 3*Me, ucm);   % naive + close correction
  U=reshape(ucm,3,Me);                                        % node-interleaved 3M -> 3 x M
  err=vecnorm(U-uex(Pe)).'/max(vecnorm(uex(Pe))); inb=dmin.'<0.1; errmax(kk)=max(err);
  fprintf('np=%2d  N=%5d/part  M=%2d  M3=%5d :  near(d<0.1) %.2e   far(d>=0.1) %.2e\n', ...
          np, Ndof, M, Me, max([err(inb);0]), max(err(~inb)));

  figure(2),clf; hold on;
  scatter3(Pe(1,:),Pe(2,:),Pe(3,:),14,log10(max(err,1e-17)),'filled');
  plot3(X{1}(1,:),X{1}(2,:),X{1}(3,:),'.','Color',[.6 .6 .6]);
  plot3(X{2}(1,:),X{2}(2,:),X{2}(3,:),'.','Color',[.6 .6 .6]);
  plot3(Y(1,:),Y(2,:),Y(3,:),'p','MarkerSize',16,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
  cb=colorbar; cb.Label.String='log_{10}|u_h-u_{exact}|/max|u_{exact}|';
  axis equal; view(35,18); grid on; colormap('jet');
  xlabel('x'); ylabel('y'); zlabel('z'); title('2-particle Stokes SLPn (S'', physmat sparse-correction): 3D velocity error (last refinement)');
  figure(2),ylim([-1 3])

end

figure(2),clf; hold on;
scatter3(Pe(1,:),Pe(2,:),Pe(3,:),14,log10(max(err,1e-17)),'filled');
plot3(X{1}(1,:),X{1}(2,:),X{1}(3,:),'.','Color',[.6 .6 .6]);
plot3(X{2}(1,:),X{2}(2,:),X{2}(3,:),'.','Color',[.6 .6 .6]);
plot3(Y(1,:),Y(2,:),Y(3,:),'p','MarkerSize',16,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
cb=colorbar; cb.Label.String='log_{10}|u_h-u_{exact}|/max|u_{exact}|';
clim([-16 -4]); axis equal; view(35,18); grid on; colormap('jet');
xlabel('x'); ylabel('y'); zlabel('z'); title('2-particle Stokes SLPn (S'', physmat sparse-correction): 3D velocity error (last refinement)');
figure(2),ylim([-1 3])

% --- convergence plot: max field error vs panel count ---
figure(1),clf; semilogy(np_vals,errmax,'o-k'); grid on; xlabel('n_p'); ylabel('max field err');
title('multi-particle Stokes SLPn (S'', physmat sparse-correction): h-refinement, p=16, M=2 n_p');

% exportgraphics(figure(1),'axissymsstok_stok_slpn_multi_physmat_convergence.png','Resolution',200)
% exportgraphics(figure(2),'axissymsstok_stok_slpn_multi_physmat_error.png','Resolution',200)

function t = applymv(x, Xall, Nall, wall, hself, hcross)
% matrix-free (-1/2 I + S') matvec on the handle framework: global stfmm3d naive traction
% (self term excluded; -1/2 jump rides in the self correction) + hself/hcross corr apply
% (each ACCUMULATES; the lab<->local R sandwich + canonical scatter happen INSIDE the module).
t = slpntracfmm(x, Xall, wall, Nall, 1e-14);
n = numel(x);
t = axpso_corr_apply_mex(hself,  n, x, n, t);                 % self correction
t = axpso_corr_apply_mex(hcross, n, x, n, t);                 % cross correction
end

function t = slpntracfmm(sig, s3dx, s3dw, NX, fmm_eps)
% naive S' traction matvec via stfmm3d at the SOURCES (self term excluded by the FMM):
% t = -p n + (grad u + grad u^T) n  (stfmm3d mu=1; S' is mu-independent).  sig INTERLOCKED.
srcinfo.sources=s3dx; srcinfo.stoklet=reshape(sig,3,[]).*s3dw;
U=stfmm3d(fmm_eps,srcinfo,3);
Gs=U.grad+permute(U.grad,[2 1 3]);                            % only the symmetric part enters
t=squeeze(sum(Gs.*reshape(NX,1,3,[]),2)) - U.pre.*NX;
t=t(:);
end
