function [At, B, F, near] = axp_offdiagphysmat_mex(iphys, ilayer, nc, nt, Xt, Nt, Ct, Rt, xct, ntc, tflag, ns, Xs, Ns, Ws, Cs, Rs, xcs, nsc, sx, snx, sws, swxp, p, np, tpan, sflag, M, iside, iclosed, mu, iforce, At, B, F, near)
iphys=double(iphys); ilayer=double(ilayer); nc=double(nc); nt=double(nt); ntc=double(ntc); tflag=double(tflag);
ns=double(ns); nsc=double(nsc); sflag=double(sflag);
p=double(p); np=double(np); M=double(M); iside=double(iside); iclosed=double(iclosed); mu=double(mu); iforce=double(iforce);
Xt=double(Xt(:)); Nt=double(Nt(:)); nt3=3*nt;                          % target 3D positions/normals (LAB), flat 3*nt
Ct=double(Ct(:)); Rt=double(Rt(:)); Cs=double(Cs(:)); Rs=double(Rs(:));
sx=double(sx(:)); snx=double(snx(:)); swxp=double(swxp(:)); sws=double(sws(:)); tpan=double(tpan(:));
if tflag ~= 0, error('axp_offdiagphysmat_mex: tflag=%d (object target) not implemented yet', tflag); end
if sflag ~= 0, error('axp_offdiagphysmat_mex: sflag=%d (curved centerline) not implemented yet', sflag); end
% ---- source panel endpoints sxlo,sxhi from the generating curve sx (dropped from the API, derived here) ----
[tg,~,~]=axt_gauss_mex(p); Le=axt_lagrange_interp_mex(p,tg,2,[-1;1]);
sxlo=complex(zeros(np,1)); sxhi=complex(zeros(np,1));
for k=1:np, ec=Le*sx((k-1)*p+(1:p)); sxlo(k)=ec(1); sxhi(k)=ec(2); end
nsp=p*np; npp1=np+1; nphi=2*M+1; ncsp=nc*nsp; Nsrc=ncsp*nphi; Np1=Nsrc+1; nnzF=nphi*nphi*ncsp; cb=ncsp*nphi; nr=nc*nt; ncns=nc*ns;
At=zeros(nr,ncns); B=complex(zeros(nr,cb)); iF=zeros(Np1,1); jF=zeros(nnzF,1); xF=complex(zeros(nnzF,1)); near=zeros(nt,1);
mex_id_ = 'axp_offdiagphysmat_r64(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i double[x], c i double[x], c i double[x], c i int64_t[x], c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i double[x], c i dcomplex[x], c i double[x], c i dcomplex[x], c i dcomplex[x], c i double[x], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io double[xx], c io dcomplex[xx], c io double[x], c io double[x], c io dcomplex[x], c io double[x])';
[At, B, iF, jF, xF, near] = AxiStokes3D_mex(mex_id_, iphys, ilayer, nc, nt, ns, Xt, Nt, Ct, Rt, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, Cs, Rs, M, iside, iclosed, mu, iforce, At, B, iF, jF, xF, near, 1, 1, 1, 1, 1, nt3, nt3, 3, 9, 1, 1, nsp, nsp, nsp, nsp, npp1, np, np, 3, 9, 1, 1, 1, 1, 1, nr, ncns, nr, cb, Np1, nnzF, nnzF, nt);
Fsq = sparse(repelem((1:Nsrc)', diff(round(iF))), jF, xF, Nsrc, Nsrc);   % source DFT kron(W,I), block-ordered [angle x comp x mer]
% ---- interlocked column reindex for the returned F (At itself comes from Fortran = real(B*F)) ----
Nmer=nsp; ncN=ncsp;
pcol=repelem((1:ns)',nc); ccx=repmat((1:nc)',ns,1);                    % interlocked col -> (surface point, component)
ii=floor((pcol-1)/Nmer)+1; jj=mod(pcol-1,Nmer)+1;                      % -> (angle i, meridian j)
bcol=(ii-1)*ncN+(ccx-1)*Nmer+jj;                                      % matching block column
F=Fsq(:,bcol);                                                        % F * sigma_interlocked = modal density
near = near > 0.5;
end

% ==== physical-space SPARSE close-correction (pass 2/3): per-panel fourier COMPUTE (Laplace SLP) ====
