function A = axss_dlp_blockmat_nmode_mex(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, mu, A, nskip, skippanidx)
if nargin<17 || isempty(nskip), nskip=0; end
if nargin<18, skippanidx=[]; end
nt=double(nt); p=double(p); np=double(np); M=double(M); iside=double(iside); iclosed=double(iclosed); mu=double(mu);
tx=double(tx(:)); sx=double(sx(:)); snx=double(snx(:)); swxp=double(swxp(:)); sxlo=double(sxlo(:)); sxhi=double(sxhi(:));
sws=double(sws(:)); tpan=double(tpan(:));
nsp=p*np; npp1=np+1; Mp1=M+1; r3t=3*nt; c3m=3*nsp*Mp1;
skippanidx=double(skippanidx(:)); skippanidx=skippanidx(1:nskip); if nskip==0, skippanidx=0; end, nsk1=max(nskip,1);
A=zeros(r3t,c3m);
mex_id_ = 'axss_dlp_blockmat_nmode_r64(c i int64_t[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i double[x], c i dcomplex[x], c i double[x], c i dcomplex[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c io dcomplex[xx], c i int64_t[x], c i int64_t[x])';
[A] = AxiStokes3D_mex(mex_id_, nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, mu, A, nskip, skippanidx, 1, nt, 1, 1, nsp, nsp, nsp, nsp, npp1, np, np, 1, 1, 1, 1, r3t, c3m, 1, nsk1);
A=reshape(A, r3t, 3*nsp, Mp1);
end

