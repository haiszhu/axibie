function A = axls_slp_blockmat_nmode_mex(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, A, nskip, skippanidx)
if nargin<16 || isempty(nskip), nskip=0; end
if nargin<17, skippanidx=[]; end
nt=double(nt); p=double(p); np=double(np); M=double(M); iside=double(iside); iclosed=double(iclosed); nskip=double(nskip);
tx=double(tx(:)); sx=double(sx(:)); snx=double(snx(:)); swxp=double(swxp(:)); sxlo=double(sxlo(:)); sxhi=double(sxhi(:));
sws=double(sws(:)); tpan=double(tpan(:));
skippanidx=double(skippanidx(:)); skippanidx=skippanidx(1:nskip); if nskip==0, skippanidx=0; end, nsk1=max(nskip,1);
nsp=p*np; npp1=np+1; Mp1=M+1; cm=nsp*Mp1;
A=zeros(nt,cm);
mex_id_ = 'axls_slp_blockmat_nmode_r64(c i int64_t[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i double[x], c i dcomplex[x], c i double[x], c i dcomplex[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c io double[xx], c i int64_t[x], c i int64_t[x])';
[A] = AxiStokes3D_mex(mex_id_, nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, A, nskip, skippanidx, 1, nt, 1, 1, nsp, nsp, nsp, nsp, npp1, np, np, 1, 1, 1, nt, cm, 1, nsk1);
A=reshape(A, nt, nsp, Mp1);
end

% SLPn 0th (target normal tnx)
