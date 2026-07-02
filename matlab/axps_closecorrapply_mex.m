function u = axps_closecorrapply_mex(nt, p, np, nang, ntcx, tcxi, idxall, S_ij, sig, iself, u)
nt=double(nt); p=double(p); np=double(np); nang=double(nang); ntcx=double(ntcx); iself=double(iself);
tcxi=double(tcxi(:)); idxall=double(idxall(:)); sig=double(sig(:)); u=double(u(:));
nsp=p*np; npp1=np+1; nsf=nsp*nang; nap=nang*p;
S_ij=double(reshape(S_ij,ntcx,nap));
if iself, nrU=nsf; else, nrU=nt; end
mex_id_ = 'axps_closecorrapply_r64(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i double[x], c i double[xx], c i double[x], c i int64_t[x], c io double[x])';
[u] = AxiStokes3D_mex(mex_id_, nt, p, np, nang, ntcx, tcxi, idxall, S_ij, sig, iself, u, 1, 1, 1, 1, 1, npp1, ntcx, ntcx, nap, nsf, 1, nrU);
end

% ============================================================
% Helpers for the pure-MATLAB n-mode block builders (utils/axls_*_blockmat_nmode.m):
% the SAME compiled gauss / lagrange-interp / carrier / sdspecialquad / split-coef routines.
% ============================================================

% Gauss-Legendre nodes/weights + spectral diff matrix on [-1,1]
