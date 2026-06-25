function A = axss_slppres_blockmat_nmode_mex(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, A)
nt=double(nt); p=double(p); np=double(np); M=double(M); iside=double(iside); iclosed=double(iclosed);
tx=double(tx(:)); sx=double(sx(:)); snx=double(snx(:)); swxp=double(swxp(:)); sxlo=double(sxlo(:)); sxhi=double(sxhi(:));
sws=double(sws(:)); tpan=double(tpan(:));
nsp=p*np; npp1=np+1; Mp1=M+1; c3m=3*nsp*Mp1;
A=zeros(nt,c3m);
mex_id_ = 'axss_slppres_blockmat_nmode_r64(c i int64_t[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i double[x], c i dcomplex[x], c i double[x], c i dcomplex[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c io dcomplex[xx])';
[A] = AxiStokes3D_mex(mex_id_, nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, A, 1, nt, 1, 1, nsp, nsp, nsp, nsp, npp1, np, np, 1, 1, 1, nt, c3m);
A=reshape(A, nt, 3*nsp, Mp1);
end

