clearvars; format short e;

% MULTI-particle Stokes SLPn (S') exterior-Neumann BVP via the Fortran operators, CLOSE-CORRECTION,
% h-refinement of np (p fixed). Mirror omega3 test_StoSLPnAxiPhysMat0: pure single-layer traction
% representation u = S[sigma], match the surface traction S'[sigma] = t_ex (the -1/2 I exterior jump
% rides in the SLPn self block; rank-deficient -> min-norm lsqminnorm). Lab node-interleaved frame:
%   SYSTEM (traction): S = naive Sto3dSLPmat-traction (FMM-replaceable) + Scorr (N self + N near diffs)
%     self  S'_kk = T_k*(Finv_k*Ab_k*F_k)*T_k'          (axp_physmat slpn)
%     near off-diag = Pi_n*(B_slpn*Fs)*T_src'           (axp_offdiagphysmat slpn, target normals)
%   EVAL (velocity): u = S[sigma] -> naive Sto3dSLPmat + SLP value close-correction (same as the SLP test).
% Solve (S+Scorr) sigma = t_ex ; evaluate the SLP velocity field.

p=16; np_vals=2:8; errmax=nan(1,numel(np_vals)); mu=1.0;
rot=@(u,th) cos(th)*eye(3)+sin(th)*[0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0]+(1-cos(th))*(u*u');
cmi=@(m) sparse([3*(1:m)-2 3*(1:m)-1 3*(1:m)].', [1:m m+(1:m) 2*m+(1:m)].', 1, 3*m, 3*m);  % comp-major(m) -> node-interleaved
Zf={@(t) 0.5*sin(t)-1i*1.0*cos(t), @(t) 1.0*sin(t)-1i*0.5*cos(t)}; ac=[0.5 1.0; 1.0 0.5];
stk=@(P,y,Fv)(1/(8*pi))*( Fv./vecnorm(P-y) + (Fv.'*(P-y)).*(P-y)./vecnorm(P-y).^3 );    % interior Stokeslet velocity
yloc={[0.10 0.05; 0 0.05; 0.20 -0.30],[0.20 -0.10; 0 0.15; 0.10 -0.10]};
Floc={[1.0 -0.5; -0.7 0.6; 0.5 -0.4],[0.6 -0.4; 0.3 -0.5; -0.7 0.5]};

for kk=1:numel(np_vals)
  np=np_vals(kk); M=2*np; nphi=2*M+1; ph=2*pi*(0:nphi-1)/nphi; iside=1; iclosed=0;
  chi_crit=cosh(log(1e13)/M);
  rng(1);

  % 2 particles, random pose; 3D node positions X{k} and 3D node normals NX{k}
  C=cell(2,1); R=cell(2,1); X=cell(2,1); NX=cell(2,1); sq=cell(2,1); geo=cell(2,1); Tc=cell(2,1); selfop=cell(2,1);
  for k=1:2
    s=[]; s.p=p; s.Z=Zf{k}; s.tpan=linspace(0,pi,np+1)'; s=quadr(s,[],'p','G');
    d=randn(3,1); C{k}=(k-1)*1.6*d/norm(d); uu=randn(3,1); uu=uu/norm(uu); R{k}=rot(uu,2*pi*rand);
    rho=real(s.x); z=imag(s.x); nr=real(s.nx); nz=imag(s.nx);
    X{k}=R{k}*[(kron(cos(ph(:)),rho)).';(kron(sin(ph(:)),rho)).';(kron(ones(nphi,1),z)).']+C{k};
    NX{k}=R{k}*[(kron(cos(ph(:)),nr)).';(kron(sin(ph(:)),nr)).';(kron(ones(nphi,1),nz)).'];   % 3D node normals
    sq{k}=axisym_to_3d_quadrature(real(s.x),imag(s.x),s.ws,nphi);
    g=[]; g.sx=s.x(:); g.snx=s.nx(:); g.sws=s.ws(:); g.swxp=s.wxp(:); g.tpan=s.tpan(:); g.sxlo=s.xlo(:); g.sxhi=s.xhi(:);
    parc=0; for kj=1:np, parc=max(parc,sum(g.sws((kj-1)*p+(1:p)))); end
    g.rnear=max(max(real(g.sx))*sqrt(2*(chi_crit-1)), 2*parc); geo{k}=g;
    [Ab,Fa,Finv,~,T]=axp_physmat_mex(2,2,3,p,np,g.sx,g.snx,g.sws,g.swxp,g.tpan,g.sxlo,g.sxhi,M,iside,iclosed,mu,R{k}); % SLPn (-1/2 I + S')
    Tc{k}=T; selfop{k}=real(T*(Finv*Ab*Fa)*T.');                 % S' lab node-interleaved self operator
  end
  Nnod=numel(geo{1}.sx)*nphi; Ndof=3*Nnod;
  Xall=[X{1} X{2}]; Nall=[NX{1} NX{2}]; wall=[sq{1}.w sq{2}.w]; Mtot=2*Nnod;
  nodeblk={1:Nnod, Nnod+1:2*Nnod}; blk={1:Ndof, Ndof+1:2*Ndof};

  % interior Stokeslets -> exact exterior velocity uex (for eval) AND exact surface traction t_ex (Neumann data)
  Y=zeros(3,0); Fg=zeros(3,0); for k=1:2, Y=[Y, R{k}*yloc{k}+C{k}]; Fg=[Fg, R{k}*Floc{k}]; end
  uex=@(P) stk(P,Y(:,1),Fg(:,1))+stk(P,Y(:,2),Fg(:,2))+stk(P,Y(:,3),Fg(:,3))+stk(P,Y(:,4),Fg(:,4));

  Picm=cmi(Mtot);
  I=[]; J=[]; V=[];                                            % COO triplets for the sparse correction operator

  % --- (2) self DIFFERENCE correction: (accurate S' self) - (naive S' self, diag 0) ---
  for k=1:2
    Pik=cmi(Nnod);
    Sns=Pik*slptrac(X{k},NX{k},X{k},sq{k}.w)*Pik.'; Sns=zeroselfblk(Sns,Nnod);
    D=selfop{k}-Sns;
    [jj,ii]=meshgrid(blk{k}, blk{k}); I=[I;ii(:)]; J=[J;jj(:)]; V=[V;D(:)];
  end

  % --- (3) off-diagonal near DIFFERENCE correction (kdtree): SLPn traction close - naive ---
  for pq=1:2
    gs=geo{pq}; other=setdiff(1:Mtot, nodeblk{pq});
    nb=~cellfun(@isempty, rangesearch(X{pq}', Xall(:,other)', gs.rnear));
    cand=other(nb); if isempty(cand), continue; end
    Pn=Xall(:,cand); phir=atan2(Pn(2,:),Pn(1,:)); txn=(hypot(Pn(1,:),Pn(2,:))+1i*Pn(3,:)).'; tphin=phir.';
    NPc=Nall(:,cand);                                          % 3D target normals -> lab-cyl (tnx meridian, tnphi azimuthal)
    nrho=NPc(1,:).*cos(phir)+NPc(2,:).*sin(phir); nphic=-NPc(1,:).*sin(phir)+NPc(2,:).*cos(phir);
    tnxc=(nrho+1i*NPc(3,:)).'; tnphic=nphic.';
    [B,Fs,nm]=axp_offdiagphysmat_mex(2,2,3,numel(cand),txn,tphin,[0;0;0],eye(3),p,np,gs.sx,gs.snx,gs.sws,gs.swxp,gs.tpan,gs.sxlo,gs.sxhi,C{pq},R{pq},M,iside,iclosed,mu,0,tnxc,tnphic);
    ncand=numel(cand); cn=find(nm); nrc=cand(cn); rr=[cn; ncand+cn; 2*ncand+cn];
    close_ni=real(cmi(numel(cn))*(B(rr,:)*Fs)*Tc{pq}.');
    naive_ni=real(cmi(numel(cn))*slptrac(Xall(:,nrc),Nall(:,nrc),X{pq},sq{pq}.w)*cmi(Nnod).');
    tdof=reshape([3*nrc-2; 3*nrc-1; 3*nrc],[],1);
    D=close_ni-naive_ni;
    [jj,ii]=meshgrid(blk{pq}, tdof); I=[I;ii(:)]; J=[J;jj(:)]; V=[V;D(:)];
  end
  Scorr=sparse(I, J, V, 2*Ndof, 2*Ndof, numel(V));

  % --- Neumann traction data + solve (pure S' rank-deficient: min-norm) ---
  trac_cm=slptrac(Xall,Nall,Y,ones(1,size(Y,2)))*[Fg(1,:).';Fg(2,:).';Fg(3,:).'];   % exact surface traction (comp-major)
  g=Picm*trac_cm;
  S = Picm*slptrac(Xall,Nall,Xall,wall)*Picm.';               % naive full S' traction, lab node-interleaved
  S = zeroselfblk(S, Mtot);
  sigma = lsqminnorm(S + Scorr, g);                           % S' rank-deficient

  % --- 3D exterior target grid (outside both particles) ---
  cmid=(C{1}+C{2})/2; gv=linspace(-3,3,22); [GX,GY,GZ]=ndgrid(gv,gv,gv); P=[GX(:)';GY(:)';GZ(:)']+cmid;
  ext=true(1,size(P,2)); dmin=inf(1,size(P,2));
  for k=1:2, loc=R{k}.'*(P-C{k}); sd=(hypot(loc(1,:),loc(2,:))/ac(k,1)).^2+(loc(3,:)/ac(k,2)).^2; ext=ext&(sd>1.05); dmin=min(dmin,sd-1); end
  Pe=P(:,ext); Me=size(Pe,2); dmin=dmin(ext);

  % --- field eval: u = S[sigma] -> naive Sto3dSLPmat + per-source near SLP-value close correction ---
  ucm=real(Sto3dSLPmat(struct('x',Pe), struct('x',Xall,'w',wall))*(Picm.'*sigma));
  for pq=1:2
    gs=geo{pq};
    cand=find(~cellfun(@isempty, rangesearch(X{pq}', Pe', gs.rnear))); if isempty(cand), continue; end
    Pn=Pe(:,cand); txn=(hypot(Pn(1,:),Pn(2,:))+1i*Pn(3,:)).'; tphin=atan2(Pn(2,:),Pn(1,:)).';
    [B,Fs,nm]=axp_offdiagphysmat_mex(2,1,3,numel(cand),txn,tphin,[0;0;0],eye(3),p,np,gs.sx,gs.snx,gs.sws,gs.swxp,gs.tpan,gs.sxlo,gs.sxhi,C{pq},R{pq},M,iside,iclosed,mu,0);
    ncand=numel(cand); cn=find(nm); sig_ni=sigma(blk{pq});
    uclose=real(B*(Fs*(Tc{pq}.'*sig_ni)));
    unaive=real(Sto3dSLPmat(struct('x',Pn), struct('x',X{pq},'w',sq{pq}.w))*(cmi(Nnod).'*sig_ni));
    for cc=1:3
      ucm((cc-1)*Me+cand(cn)) = ucm((cc-1)*Me+cand(cn)) - unaive((cc-1)*ncand+cn) + uclose((cc-1)*ncand+cn);
    end
  end
  U=[ucm(1:Me).'; ucm(Me+1:2*Me).'; ucm(2*Me+1:3*Me).'];
  err=vecnorm(U-uex(Pe)).'/max(vecnorm(uex(Pe))); inb=dmin.'<0.1; errmax(kk)=max(err);
  fprintf('np=%2d  N=%5d/part  M=%2d  M3=%5d :  near(d<0.1) %.2e   far(d>=0.1) %.2e\n', ...
          np, Ndof, M, Me, max([err(inb);0]), max(err(~inb)));

  figure(2),clf; hold on;
  scatter3(Pe(1,:),Pe(2,:),Pe(3,:),14,log10(max(err,1e-17)),'filled');
  plot3(X{1}(1,:),X{1}(2,:),X{1}(3,:),'.','Color',[.6 .6 .6]);
  plot3(X{2}(1,:),X{2}(2,:),X{2}(3,:),'.','Color',[.6 .6 .6]);
  plot3(Y(1,:),Y(2,:),Y(3,:),'p','MarkerSize',16,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
  cb=colorbar; cb.Label.String='log_{10}|u_h-u_{exact}|/max|u_{exact}|';
  clim([-16 -4]); axis equal; view(35,18); grid on; colormap('jet');
  xlabel('x'); ylabel('y'); zlabel('z'); title('2-particle Stokes SLPn (S'', close-correction): 3D velocity error (last refinement)');
end

% --- convergence plot: max field error vs panel count ---
figure(1),clf; semilogy(np_vals,errmax,'o-k'); grid on; xlabel('n_p'); ylabel('max field err');
title('multi-particle Stokes SLPn (S'', close-correction): h-refinement, p=16, M=2 n_p');

% cmp = getPyPlot_cMap('rainbow', [], [], '"/Users/hzhu/.pyenv/versions/3.11.13/bin/python"');
% colormap(cmp)
% 
% exportgraphics(figure(1),'axissymsstok_stok_slpn_multi_convergence.png','Resolution',200)
% exportgraphics(figure(2),'axissymsstok_stok_slpn_multi_error.png','Resolution',200)

function T = slptrac(Xt, Nt, xs, ws)
% S' (single-layer traction) operator: target traction matrix from a source struct.  Needs target normals Nt.
[~,T] = Sto3dSLPmat(struct('x',Xt,'nx',Nt), struct('x',xs,'w',ws)); T = real(T);
end

function A = zeroselfblk(A, m)
% zero the m node-diagonal 3x3 blocks of a 3m x 3m node-interleaved matrix (singular self-interaction)
N3 = 3*m; ri = (1:3).' + 3*(0:m-1);
lin = (reshape(ri,3,1,m)-1)*N3 + reshape(ri,1,3,m);
A(lin(:)) = 0;
end
