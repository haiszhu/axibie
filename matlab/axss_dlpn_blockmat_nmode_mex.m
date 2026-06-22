function A = axss_dlpn_blockmat_nmode_mex(nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, mu, A)
nt=double(nt); p=double(p); np=double(np); M=double(M); iside=double(iside); iclosed=double(iclosed); mu=double(mu);
tx=double(tx(:)); tnx=double(tnx(:)); sx=double(sx(:)); snx=double(snx(:)); swxp=double(swxp(:)); sxlo=double(sxlo(:)); sxhi=double(sxhi(:));
sws=double(sws(:)); tpan=double(tpan(:));
nsp=p*np; npp1=np+1; Mp1=M+1; r3t=3*nt; c3m=3*nsp*Mp1;
A=zeros(r3t,c3m);
mex_id_ = 'axss_dlpn_blockmat_nmode_r64(c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i double[x], c i dcomplex[x], c i double[x], c i dcomplex[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c io dcomplex[xx])';
[A] = AxiStokes3D_mex(mex_id_, nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, mu, A, 1, nt, nt, 1, 1, nsp, nsp, nsp, nsp, npp1, np, np, 1, 1, 1, 1, r3t, c3m);
A=reshape(A, r3t, 3*nsp, Mp1);
end

% DLP n-mode
