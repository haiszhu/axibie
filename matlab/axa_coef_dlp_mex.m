function [C1,C2,C3,C4,C5] = axa_coef_dlp_mex(nt, nq, M, mu, zr, zi, yr, yi, nvr, nvi, C1, C2, C3, C4, C5)
nt=double(nt); nq=double(nq); M=double(M); mu=double(mu);
zr=double(zr(:)); zi=double(zi(:)); yr=double(yr(:)); yi=double(yi(:)); nvr=double(nvr(:)); nvi=double(nvi(:));
Mp1=M+1; NN=9*nt*nq;
C1=zeros(Mp1,NN); C2=zeros(Mp1,NN); C3=zeros(Mp1,NN); C4=zeros(Mp1,NN); C5=zeros(Mp1,NN);
mex_id_ = 'axa_coef_dlp_r64(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c io double[xx], c io double[xx], c io dcomplex[xx], c io dcomplex[xx], c io dcomplex[xx])';
[C1, C2, C3, C4, C5] = AxiStokes3D_mex(mex_id_, nt, nq, M, mu, zr, zi, yr, yi, nvr, nvi, C1, C2, C3, C4, C5, 1, 1, 1, 1, nt, nt, nq, nq, nq, nq, Mp1, NN, Mp1, NN, Mp1, NN, Mp1, NN, Mp1, NN);
end

% SLPn split coefficients (MER + SWIRL).  TARGET normal nxr,nxi.  MER: C1,C2 real, C3,C4 cplx.
% SWIRL (unit n_theta): C1s,C2s real, C3s cplx (C4s=0).  Each (M+1)x(9 nt nq); reshape (M+1,3nt,3nq).
