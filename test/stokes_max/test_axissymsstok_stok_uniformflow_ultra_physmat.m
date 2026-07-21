% /Applications/ParaView-6.1.0.app/Contents/bin/pvpython pv_flow_recipe.py \
%     --stem axissymsstok_uniformflow_ultra --k 27  --seed-radius-frac 0.9 --slice-opacity 0.3 --seed-points 200
% Linux (ParaView 5.11.2 on PATH):
% /usr/bin/pvpython pv_flow_recipe.py \
%     --stem axissymsstok_uniformflow_ultra --k 27  --seed-radius-frac 0.9 --slice-opacity 0.3 --seed-points 200
clearvars; format short e;
addpath('../../../axibie/utils');                            % axibie sibling repo (Github/axibie next to Github/AxiStokes3D)
addpath('../../utils');
addpath('../../matlab');
addpath('../../external/fmm3d/matlab');                       % stfmm3d

% Uniform BACKGROUND FLOW past a Kside^3 lattice of FIXED rigid spheroids.  GEOMETRY SETUP copied
% from test_axissymsstok_stok_slpn_bvp_ultra_physmat2.m (alternating prolate/oblate shapes, random
% rotations, lattice spacing Lsp, gate=2.0, Rbound+rnear near-zone).  The ONLY change vs the SLPn
% benchmark is the physics: a Dirichlet (velocity) drive, so a velocity combined-field (D+S)[sigma]
% is used (the ilayer=3 DLP + ilayer=1 SLP solver-module machinery of
% test_axissymsstok_stok_dlp_bvp_max_physmat.m), and the RHS is the uniform ambient Uinf.
%
% No-slip on each FIXED body -> total velocity 0 on the surface (exterior jump +1/2 I rides in the
% DLP self correction; S regularizes the DLP rigid-body null space -> SECOND-KIND full-rank, gmres
% needs NO deflation / NO preconditioner):
%     (1/2 I + D + S)[sigma] = -Uinf
% Total lab-frame field  u(x) = Uinf + (D+S)[sigma](x)   ( -> 0 on the surfaces, -> Uinf at infinity).
% Four corr_setup handles {D,S} x {self,cross} applied SEPARATELY (no combine bookkeeping).  Naive
% (D+S) velocity is ONE stfmm3d call (stoklet + strslet/strsvec, pot).  NO exact solution / NO
% convergence: output = the ndgrid velocity field (.vti) + particle surfaces (.vtp) for ParaView.
% Drag on body k (net Stokeslet strength) = sum sigma_k .* w_k.  mu=1, node-interleaved lab frame.

p=16; np=4; mu=1.0; gate=2.0; iside=1; iclosed=0;
M=2*np; nphi=2*M+1;
Kside=3; K=Kside^3;                                           % ultra_physmat2 uses Kside=15 (K=3375, 512GB); 3 (K=27) runs/visualizes locally
% Kside=8; K=Kside^3;
p=8; np=4; Kside=10; K=Kside^3;
fmmeps=1e-12;
Uinf=[1;0;0]; shr=0.5;                                        % uniform velocity Uinf / simple-shear rate shr (settable)
% --- ambient background flow  ubg(P): 3 x n, MUST be a Stokes solution; used in BOTH the RHS and the eval ---
ubg = @(P) repmat(Uinf,1,size(P,2));                        % UNIFORM               (u -> Uinf at infinity)
% ubg = @(P) [shr*P(3,:); zeros(2,size(P,2))];                % SHEAR, gradient along z (screen-horizontal, depth)
% ubg = @(P) [zeros(2,size(P,2)); shr*P(1,:)];                % SHEAR, gradient along x (screen-horizontal)
% ubg = @(P) [shr*P(2,:); zeros(2,size(P,2))];                  % SHEAR, gradient along y (screen UP-DOWN): flow x, u_x = shr*y  <- active
rot=@(u,th) cos(th)*eye(3)+sin(th)*[0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0]+(1-cos(th))*(u*u');
Zf={@(t) 0.5*sin(t)-1i*1.0*cos(t), @(t) 1.0*sin(t)-1i*0.5*cos(t)};   % prolate / oblate
ac=[0.5 1.0; 1.0 0.5]; Rbound=1.0;

% ---- per-shape meridian geometry (self passes are PER PARTICLE) ----  [ultra copy]
geo=cell(2,1); sq=cell(2,1);
for sh=1:2
  s=[]; s.p=p; s.Z=Zf{sh}; s.tpan=linspace(0,pi,np+1)'; s=quadr(s,[],'p','G');
  g=[]; g.sx=s.x(:); g.snx=s.nx(:); g.sws=s.ws(:); g.swxp=s.wxp(:); g.tpan=s.tpan(:);
  g.rho=real(s.x); g.z=imag(s.x); g.nr=real(s.nx); g.nz=imag(s.nx); geo{sh}=g;
  sq{sh}=axisym_to_3d_quadrature(g.rho,g.z,s.ws,g.nr,g.nz,nphi);
end
N=numel(geo{1}.sx); Nnod=N*nphi;
rnear=gate*sqrt(max(arrayfun(@(j) sum(geo{1}.sws((j-1)*p+(1:p))), 1:np)));   % near-zone width
Lsp=2*Rbound + 0.4;                                          % lattice spacing (surface gap 0.4, ultra convention)
fprintf('K=%d  p=%d np=%d nphi=%d (Nnod=%d/particle, dof=%d)  gate=%.1f  spacing=%.2f\n', ...
        K, p, np, nphi, Nnod, 3*K*Nnod, gate, Lsp);
rng(1);

% 1. lattice: alternating shapes, random rotations  [ultra copy]
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

% 2. uniform background-flow BC: no-slip on the FIXED bodies -> (1/2 I + D + S)[sigma] = -Uinf
g_di = reshape(-ubg(Xall), [], 1);                          % no-slip: (1/2 I + D + S)[sigma] = -ubg at the surface nodes

% flat per-particle meridian geometry (consumed by all corr_setup calls)  [ultra copy]
sxs=zeros(N,K); snxs=sxs; swss=zeros(N,K); swxps=zeros(N,K); tpans=zeros(np+1,K);
for k=1:K
  gs=geo{ty(k)};
  sxs(:,k)=gs.sx; snxs(:,k)=gs.snx; swss(:,k)=gs.sws; swxps(:,k)=gs.swxp; tpans(:,k)=gs.tpan;
end
pv=p*ones(K,1); npv=np*ones(K,1); pmv=M*ones(K,1);
geomoff=1+(0:K)'*N; tpanoff=1+(0:K)'*(np+1); targoff=1+(0:K)'*Nnod;
nsx=K*N; ntpan=K*(np+1);
sxf=sxs(:); snxf=snxs(:); swsf=swss(:); swxpf=swxps(:); tpanf=tpans(:);

% 3. corrections: FOUR corr_setup handles {D (ilayer=3, stresslet + 1/2 jump), S (ilayer=1, Stokeslet)}
%    x {self (iinter=1), cross (iinter=2)} -- applied SEPARATELY in the matvec (no combine).
hDself  = axpso_corr_setup_mex(2, 3, 1, 0, mu, 1, K, pv, npv, pmv, iside, iclosed, gate, ...
            geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
            Rbound+rnear, K*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);
hSself  = axpso_corr_setup_mex(2, 1, 1, 0, mu, 1, K, pv, npv, pmv, iside, iclosed, gate, ...
            geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
            Rbound+rnear, K*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);
hDcross = axpso_corr_setup_mex(2, 3, 1, 0, mu, 2, K, pv, npv, pmv, iside, iclosed, gate, ...
            geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
            Rbound+rnear, K*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);
hScross = axpso_corr_setup_mex(2, 1, 1, 0, mu, 2, K, pv, npv, pmv, iside, iclosed, gate, ...
            geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
            Rbound+rnear, K*Nnod, targoff, Xall, Nall, cat(3,R{:}), [C{:}], 0);

% 4. matrix-free (1/2 I + D + S) matvec; SECOND-KIND full-rank -> no deflation, no preconditioner
applyA=@(x) applyDS(x, fmmeps, Xall, Nall, wall, hDself, hSself, hDcross, hScross);

% 5. solve
% [sigma,flag,relres,iter]=gmres(applyA, g_di, [], 1e-9, 200);
[sigma,flag,relres,iter]=gmres(applyA, g_di, [], 1e-4, 200);
fprintf('gmres: flag=%d  iters=%d  relres=%.3e\n', flag, iter(end), relres);
sig3=reshape(sigma,3,K*Nnod); Fk=zeros(3,K);                % drag body k = net Stokeslet strength = sum sigma_k .* w_k
for k=1:K, idx=(k-1)*Nnod+(1:Nnod); Fk(:,k)=sum(sig3(:,idx).*wall(idx),2); end
fprintf('drag: total F=[%.3f %.3f %.3f], per-body |F| in [%.4f, %.4f]\n', sum(Fk,2), min(vecnorm(Fk)), max(vecnorm(Fk)));

% 6. off-surface ndgrid (x-fastest = VTK order); exact spheroid inside-filter  [physmat5 grid]
ng=48; margin=0.9;
ng=96; margin=4;      % margin = clear-fluid halo (abs units, ~Rbound) around the particle bed;
                      % ~2 Rbound so approach flow / wake are visible, not buried in the cluster
ng=144; margin=4;                      
glo=-Lsp*(Kside-1)/2 - Rbound - margin;
xs=linspace(glo,-glo,ng);
[xg,yg,zg]=ndgrid(xs,xs,xs); Pg=[xg(:).'; yg(:).'; zg(:).'];
inside=false(1,size(Pg,2)); dall=inf(1,size(Pg,2));
for k=1:K
  loc=R{k}.'*(Pg-C{k}); sd=(hypot(loc(1,:),loc(2,:))/ac(ty(k),1)).^2+(loc(3,:)/ac(ty(k),2)).^2;
  inside=inside | (sd<1.0); dall=min(dall,sqrt(sd)-1);
end
ext=~inside; Pe=Pg(:,ext); Me=size(Pe,2);
fprintf('grid: %d^3, %d exterior targets (h=%.3f)\n', ng, Me, xs(2)-xs(1));

% 7. field eval: naive (D+S) velocity = ONE stfmm3d pottarg (stoklet + strslet/strsvec) + {D,S}
%    eval-handle applies (targoff=ones -> no own-block drop); total u = Uinf + disturbance.
sig3w=reshape(sigma,3,[]).*wall;
srcinfo=[]; srcinfo.sources=Xall; srcinfo.stoklet=sig3w; srcinfo.strslet=sig3w; srcinfo.strsvec=Nall;
Uf=stfmm3d(fmmeps,srcinfo,0,Pe,1);
uc=Uf.pottarg(:);
hDeval = axpso_corr_setup_mex(2, 3, 1, 0, mu, 2, K, pv, npv, pmv, iside, iclosed, gate, ...
           geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
           Rbound+rnear, Me, ones(K+1,1), Pe, Pe, cat(3,R{:}), [C{:}], 0);
hSeval = axpso_corr_setup_mex(2, 1, 1, 0, mu, 2, K, pv, npv, pmv, iside, iclosed, gate, ...
           geomoff, tpanoff, nsx, ntpan, sxf, snxf, swsf, swxpf, tpanf, ...
           Rbound+rnear, Me, ones(K+1,1), Pe, Pe, cat(3,R{:}), [C{:}], 0);
uc=axpso_corr_apply_mex(hDeval, 3*K*Nnod, sigma, 3*Me, uc);
uc=axpso_corr_apply_mex(hSeval, 3*K*Nnod, sigma, 3*Me, uc);
Uh=zeros(3,size(Pg,2)); Uh(:,ext)=ubg(Pe)+reshape(uc,3,[]); % total lab-frame velocity = ambient + disturbance (0 inside)
far=ext & (dall>2);                                         % consistency (no exact soln): disturbance -> 0, i.e. u -> ubg far
fprintf('disturbance decay (d>2, %d pts): max|u - ubg| = %.2e\n', sum(far), max(vecnorm(Uh(:,far)-ubg(Pg(:,far)))));

% 8. save the velocity grid (.vti) + particle surfaces (.vtp) for ParaView (arrays: u, umag, inside)  [physmat5]
pdv=[]; pdv.u=Uh; pdv.umag=vecnorm(Uh); pdv.inside=double(inside);
write_vti_grid(sprintf('axissymsstok_uniformflow_ultra_K%d.vti',K), ...
               [glo glo glo], (xs(2)-xs(1))*[1 1 1], [ng ng ng], pdv);
nmer=48; naz=64; tm=linspace(0,pi,nmer); phd=2*pi*(0:naz-1)/naz;
VV=zeros(3,K*nmer*naz); shp=zeros(1,K*nmer*naz);
F0=zeros((nmer-1)*naz,4); c=0;
for aa=1:naz
  a2=mod(aa,naz)+1;                                          % azimuthal wrap
  for j=1:nmer-1, c=c+1; F0(c,:)=[(aa-1)*nmer+j, (aa-1)*nmer+j+1, (a2-1)*nmer+j+1, (a2-1)*nmer+j]; end
end
for k=1:K
  zc=Zf{ty(k)}(tm); rho=real(zc); zz=imag(zc);
  Xk=[kron(cos(phd(:)),rho(:)).'; kron(sin(phd(:)),rho(:)).'; kron(ones(naz,1),zz(:)).'];
  idxv=(k-1)*nmer*naz+(1:nmer*naz); VV(:,idxv)=R{k}*Xk+C{k}; shp(idxv)=ty(k);
end
FF=repmat(F0,K,1)+repelem((0:K-1)'*(nmer*naz), size(F0,1), 1)*ones(1,4);
pds=[]; pds.shape=shp;
write_vtp_mesh(sprintf('axissymsstok_uniformflow_ultra_K%d.vtp',K), VV, FF, pds);
fprintf('saved axissymsstok_uniformflow_ultra_K%d.{vti,vtp}\n', K);

% -------------------------------------------------------------------------------------------
function y=applyDS(x, fmmeps, Xall, Nall, wall, hDself, hSself, hDcross, hScross)
% matrix-free (+1/2 I + D + S) matvec: ONE stfmm3d call (stoklet + strslet/strsvec) = naive (D+S)
% velocity at the sources (self i~=j excluded; +1/2 jump rides in the DLP self correction) + the
% four handle corr applies (ACCUMULATE, R sandwich INSIDE the module).  x node-interleaved.
sig3=reshape(x,3,[]).*wall;
srcinfo.sources=Xall; srcinfo.stoklet=sig3; srcinfo.strslet=sig3; srcinfo.strsvec=Nall;
U=stfmm3d(fmmeps,srcinfo,1);
y=U.pot(:);
n=numel(x);
y=axpso_corr_apply_mex(hDself,  n, x, n, y);
y=axpso_corr_apply_mex(hSself,  n, x, n, y);
y=axpso_corr_apply_mex(hDcross, n, x, n, y);
y=axpso_corr_apply_mex(hScross, n, x, n, y);
end
