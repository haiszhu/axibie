function [A, Ab, F, Finv, Abinv, T] = axp_physmat_mex(iphys, ilayer, nc, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, mu, R)
iphys=double(iphys); ilayer=double(ilayer); nc=double(nc); p=double(p); np=double(np); M=double(M); iside=double(iside); iclosed=double(iclosed); mu=double(mu);
sx=double(sx(:)); snx=double(snx(:)); swxp=double(swxp(:)); sxlo=double(sxlo(:)); sxhi=double(sxhi(:));
sws=double(sws(:)); tpan=double(tpan(:)); R=double(R(:));
nsp=p*np; npp1=np+1; nphi=2*M+1; r=nc*nsp; Nfull=r*nphi; Np1=Nfull+1;
nnzA=nphi*r*r; nnzF=nphi*nphi*r; nnzT=nc*nc*nsp*nphi;
A=zeros(Nfull,Nfull);                                                 % assembled SELF operator, LAB node-interleaved (real)
iAb=zeros(Np1,1); jAb=zeros(nnzA,1); xAb=complex(zeros(nnzA,1));
iAi=zeros(Np1,1); jAi=zeros(nnzA,1); xAi=complex(zeros(nnzA,1));
iF=zeros(Np1,1); jF=zeros(nnzF,1); xF=complex(zeros(nnzF,1));
iFi=zeros(Np1,1); jFi=zeros(nnzF,1); xFi=complex(zeros(nnzF,1));
iT=zeros(Np1,1); jT=zeros(nnzT,1); xT=zeros(nnzT,1);
mex_id_ = 'axp_physmat_r64(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i double[x], c i dcomplex[x], c i double[x], c i dcomplex[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i double[x], c io double[xx], c io double[x], c io double[x], c io dcomplex[x], c io double[x], c io double[x], c io dcomplex[x], c io double[x], c io double[x], c io dcomplex[x], c io double[x], c io double[x], c io dcomplex[x], c io double[x], c io double[x], c io double[x])';
[A, iAb, jAb, xAb, iAi, jAi, xAi, iF, jF, xF, iFi, jFi, xFi, iT, jT, xT] = AxiStokes3D_mex(mex_id_, iphys, ilayer, nc, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, mu, R, A, iAb, jAb, xAb, iAi, jAi, xAi, iF, jF, xF, iFi, jFi, xFi, iT, jT, xT, 1, 1, 1, 1, 1, nsp, nsp, nsp, nsp, npp1, np, np, 1, 1, 1, 1, 9, Nfull, Nfull, Np1, nnzA, nnzA, Np1, nnzA, nnzA, Np1, nnzF, nnzF, Np1, nnzF, nnzF, Np1, nnzT, nnzT);
% CSR (ia=row_ptr, ja=col, a=val) -> MATLAB sparse (triplet): expand row_ptr to row indices
rows=(1:Nfull)';
Ab    = sparse(repelem(rows, diff(round(iAb))), jAb, xAb, Nfull, Nfull);
Abinv = sparse(repelem(rows, diff(round(iAi))), jAi, xAi, Nfull, Nfull);
F     = sparse(repelem(rows, diff(round(iF))),  jF,  xF,  Nfull, Nfull);
Finv  = sparse(repelem(rows, diff(round(iFi))), jFi, xFi, Nfull, Nfull);
T     = sparse(repelem(rows, diff(round(iT))),  jT,  xT,  Nfull, Nfull);
end

% ==== physical-space OFF-DIAGONAL operator (AxiOffDiagPhysMat): one source -> arbitrary targets ====
