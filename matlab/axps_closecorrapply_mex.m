function u = axps_closecorrapply_mex(nc, nt, p, np, nang, ntcx, tcxi, idxall, S_ij, sig, iself, u)
nc=double(nc); nt=double(nt); p=double(p); np=double(np); nang=double(nang); ntcx=double(ntcx); iself=double(iself);
tcxi=double(tcxi(:)); idxall=double(idxall(:)); sig=double(sig(:)); u=double(u(:));
nsp=p*np; npp1=np+1; ncnsf=nc*nsp*nang; ncnr=nc*ntcx; ncnap=nc*nang*p;
S_ij=double(reshape(S_ij,ncnr,ncnap));
if iself, nrU=ncnsf; else, nrU=nc*nt; end
mex_id_ = 'axps_closecorrapply_r64(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i double[x], c i double[xx], c i double[x], c i int64_t[x], c io double[x])';
[u] = AxiStokes3D_mex(mex_id_, nc, nt, p, np, nang, ntcx, tcxi, idxall, S_ij, sig, iself, u, 1, 1, 1, 1, 1, 1, npp1, ntcx, ncnr, ncnap, ncnsf, 1, nrU);
end

% ==== physical-space SPARSE close-correction (pass 2/3): Stokes SLP velocity COMPUTE (nc=3) ====
