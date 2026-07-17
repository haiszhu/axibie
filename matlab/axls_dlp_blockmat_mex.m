function A = axls_dlp_blockmat_mex(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, A)
nt=double(nt); p=double(p); np=double(np); iside=double(iside); iclosed=double(iclosed);
tx=double(tx(:)); sx=double(sx(:)); snx=double(snx(:)); swxp=double(swxp(:)); sxlo=double(sxlo(:)); sxhi=double(sxhi(:));
sws=double(sws(:)); tpan=double(tpan(:));
nsp=p*np; npp1=np+1;
A=zeros(nt,nsp);
mex_id_ = 'axls_dlp_blockmat_r64(c i int64_t[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i double[x], c i dcomplex[x], c i double[x], c i dcomplex[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c io double[xx])';
[A] = AxiStokes3D_mex(mex_id_, nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, A, 1, nt, 1, 1, nsp, nsp, nsp, nsp, npp1, np, np, 1, 1, nt, nsp);
end

% DLP n-mode (optional trailing nskip, skippanidx: coarse panels 1..np whose p columns stay zero)
