function [ntcxs, nbytes] = axpso_stokslpn_selfsetup_mex(K, p, np, nang, pmodes, iside, iclosed, gate, mu, sx, snx, sws, swxp, tpan, ntcxs, nbytes)
K=double(K); p=double(p); np=double(np); nang=double(nang); pmodes=double(pmodes);
iside=double(iside); iclosed=double(iclosed); gate=double(gate); mu=double(mu);
nsp=p*np; npp1=np+1;
sx=reshape(double(sx),nsp,K); snx=reshape(double(snx),nsp,K); swxp=reshape(double(swxp),nsp,K);
sws=reshape(double(sws),nsp,K); tpan=reshape(double(tpan),npp1,K);
ntcxs=double(ntcxs(:)); nbytes=double(nbytes);
mex_id_ = 'axpso_stokslpn_selfsetup_r64(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i double[x], c i dcomplex[xx], c i dcomplex[xx], c i double[xx], c i dcomplex[xx], c i double[xx], c io double[x], c io double[x])';
[ntcxs, nbytes] = AxiStokes3D_mex(mex_id_, K, p, np, nang, pmodes, iside, iclosed, gate, mu, sx, snx, sws, swxp, tpan, ntcxs, nbytes, 1, 1, 1, 1, 1, 1, 1, 1, 1, nsp, K, nsp, K, nsp, K, nsp, K, npp1, K, K, 1);
end

