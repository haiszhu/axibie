function [ntcxx, nuo, npair, nbytes] = axpso_stokslpn_crosssetup_mex(K, p, np, nang, pmodes, iside, iclosed, gate, mu, rball, Xall, Nall, Rm, Cc, sx, snx, sws, swxp, tpan, ntcxx, nuo, npair, nbytes)
K=double(K); p=double(p); np=double(np); nang=double(nang); pmodes=double(pmodes);
iside=double(iside); iclosed=double(iclosed); gate=double(gate); mu=double(mu); rball=double(rball);
nsp=p*np; npp1=np+1; ktot=nsp*nang*K; k9=9*K; k3=3*K;
Xall=double(reshape(Xall,3,ktot)); Nall=double(reshape(Nall,3,ktot));
Rm=double(reshape(Rm,9,K)); Cc=double(reshape(Cc,3,K));
sx=reshape(double(sx),nsp,K); snx=reshape(double(snx),nsp,K); swxp=reshape(double(swxp),nsp,K);
sws=reshape(double(sws),nsp,K); tpan=reshape(double(tpan),npp1,K);
ntcxx=double(ntcxx(:)); nuo=double(nuo(:)); npair=double(npair); nbytes=double(nbytes);
mex_id_ = 'axpso_stokslpn_crosssetup_r64(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i double[x], c i double[x], c i double[xx], c i double[xx], c i double[x], c i double[x], c i dcomplex[xx], c i dcomplex[xx], c i double[xx], c i dcomplex[xx], c i double[xx], c io double[x], c io double[x], c io double[x], c io double[x])';
[ntcxx, nuo, npair, nbytes] = AxiStokes3D_mex(mex_id_, K, p, np, nang, pmodes, iside, iclosed, gate, mu, rball, Xall, Nall, Rm, Cc, sx, snx, sws, swxp, tpan, ntcxx, nuo, npair, nbytes, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, ktot, 3, ktot, k9, k3, nsp, K, nsp, K, nsp, K, nsp, K, npp1, K, K, K, 1, 1);
end

