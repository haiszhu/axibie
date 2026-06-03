function [C1,C2,C3,C4,C1s,C2s,C3s] = axa_coef_slpn_mex(nt, nq, M, mu, zr, zi, yr, yi, nxr, nxi, C1, C2, C3, C4, C1s, C2s, C3s)
nt=double(nt); nq=double(nq); M=double(M); mu=double(mu);
zr=double(zr(:)); zi=double(zi(:)); yr=double(yr(:)); yi=double(yi(:)); nxr=double(nxr(:)); nxi=double(nxi(:));
Mp1=M+1; NN=9*nt*nq;
C1=zeros(Mp1,NN); C2=zeros(Mp1,NN); C3=zeros(Mp1,NN); C4=zeros(Mp1,NN); C1s=zeros(Mp1,NN); C2s=zeros(Mp1,NN); C3s=zeros(Mp1,NN);
mex_id_ = 'axa_coef_slpn_r64(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c io double[xx], c io double[xx], c io dcomplex[xx], c io dcomplex[xx], c io double[xx], c io double[xx], c io dcomplex[xx])';
[C1, C2, C3, C4, C1s, C2s, C3s] = AxiStokes3D_mex(mex_id_, nt, nq, M, mu, zr, zi, yr, yi, nxr, nxi, C1, C2, C3, C4, C1s, C2s, C3s, 1, 1, 1, 1, nt, nt, nq, nq, nt, nt, Mp1, NN, Mp1, NN, Mp1, NN, Mp1, NN, Mp1, NN, Mp1, NN, Mp1, NN);
end

% DLP close-eval operator block (port of StokDLPAxiModalSpecialquad).  Complex inputs split
% real/imag; host supplies upsampled weights gw (q=2p) and the q->p Legendre projection P (q x p).
% W is (M+1)x(9 nt p); reshape to (M+1,3nt,3p) on return.
