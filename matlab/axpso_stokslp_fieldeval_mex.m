function u = axpso_stokslp_fieldeval_mex(K, Me, p, np, nang, pmodes, iside, iclosed, gate, mu, rball, Pe, sigma, Rm, Cc, sx, snx, sws, swxp, tpan, u)
K=double(K); Me=double(Me); p=double(p); np=double(np); nang=double(nang); pmodes=double(pmodes);
iside=double(iside); iclosed=double(iclosed); gate=double(gate); mu=double(mu); rball=double(rball);
nsp=p*np; npp1=np+1; me3=3*Me; nsig=3*K*nsp*nang; k9=9*K; k3=3*K;
Pe=double(reshape(Pe,3,Me)); sigma=double(sigma(:)); Rm=double(reshape(Rm,9,K)); Cc=double(reshape(Cc,3,K));
sx=reshape(double(sx),nsp,K); snx=reshape(double(snx),nsp,K); swxp=reshape(double(swxp),nsp,K);
sws=reshape(double(sws),nsp,K); tpan=reshape(double(tpan),npp1,K); u=double(u(:));
mex_id_ = 'axpso_stokslp_fieldeval_r64(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i double[x], c i double[x], c i double[xx], c i double[x], c i double[x], c i double[x], c i dcomplex[xx], c i dcomplex[xx], c i double[xx], c i dcomplex[xx], c i double[xx], c io double[x])';
[u] = AxiStokes3D_mex(mex_id_, K, Me, p, np, nang, pmodes, iside, iclosed, gate, mu, rball, Pe, sigma, Rm, Cc, sx, snx, sws, swxp, tpan, u, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, Me, nsig, k9, k3, nsp, K, nsp, K, nsp, K, nsp, K, npp1, K, me3);
end

