function A = axls_dlpn_blockmat_nmode_mex(nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, A)
nt=double(nt); p=double(p); np=double(np); M=double(M); iside=double(iside); iclosed=double(iclosed);
tx=double(tx(:)); tnx=double(tnx(:)); sx=double(sx(:)); snx=double(snx(:)); swxp=double(swxp(:)); sxlo=double(sxlo(:)); sxhi=double(sxhi(:));
sws=double(sws(:)); tpan=double(tpan(:));
nsp=p*np; npp1=np+1; Mp1=M+1; cm=nsp*Mp1;
A=zeros(nt,cm);
mex_id_ = 'axls_dlpn_blockmat_nmode_r64(c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i double[x], c i dcomplex[x], c i double[x], c i dcomplex[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c io double[xx])';
[A] = AxiStokes3D_mex(mex_id_, nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, A, 1, nt, nt, 1, 1, nsp, nsp, nsp, nsp, npp1, np, np, 1, 1, 1, nt, cm);
A=reshape(A, nt, nsp, Mp1);
end
