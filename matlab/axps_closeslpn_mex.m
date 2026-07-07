function [S_ij, idxall] = axps_closeslpn_mex(nt, tx, t3dx, tn3dx, p, np, nang, sx, snx, sws, swxp, tpan, gate, s3dx, s3dnx, s3dw, pmodes, iside, iclosed, ntcx, tcxi, S_ij, idxall)
nt=double(nt); p=double(p); np=double(np); nang=double(nang); pmodes=double(pmodes);
iside=double(iside); iclosed=double(iclosed); ntcx=double(ntcx); gate=double(gate);
tx=double(tx(:)); sx=double(sx(:)); snx=double(snx(:)); swxp=double(swxp(:));
sws=double(sws(:)); tpan=double(tpan(:)); tcxi=double(tcxi(:)); s3dw=double(s3dw(:));
nsp=p*np; npp1=np+1; nsf=nsp*nang; nt3=3*nt; ns3=3*nsf; nap=nang*p;
t3dx=double(reshape(t3dx,3,nt)); tn3dx=double(reshape(tn3dx,3,nt)); s3dx=double(reshape(s3dx,3,nsf)); s3dnx=double(reshape(s3dnx,3,nsf));
S_ij=double(reshape(S_ij,ntcx,nap)); idxall=double(idxall(:));
mex_id_ = 'axps_closeslpn_r64(c i int64_t[x], c i dcomplex[x], c i double[x], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i double[x], c i dcomplex[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c io double[xx], c io double[x])';
[S_ij, idxall] = AxiStokes3D_mex(mex_id_, nt, tx, t3dx, tn3dx, p, np, nang, sx, snx, sws, swxp, tpan, gate, s3dx, s3dnx, s3dw, pmodes, iside, iclosed, ntcx, tcxi, S_ij, idxall, 1, nt, nt3, nt3, 1, 1, 1, nsp, nsp, nsp, nsp, npp1, 1, ns3, ns3, nsf, 1, 1, 1, 1, npp1, ntcx, nap, ntcx);
end

% ==== physical-space SPARSE close-correction (pass 2/3): per-panel fourier COMPUTE (Laplace DLP) ====
