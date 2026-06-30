function [C1, C2, C3] = axsk_slppres_coef_nmode_mex(nt, tx, nq, sx, M)
nt=double(nt); nq=double(nq); M=double(M); tx=double(tx(:)); sx=double(sx(:));
Mp1=M+1; c3m=3*nq*Mp1;
C1=complex(zeros(nt,c3m)); C2=complex(zeros(nt,c3m)); C3=complex(zeros(nt,c3m));
mex_id_ = 'axsk_slppres_coef_nmode_r64(c i int64_t[x], c i dcomplex[x], c i int64_t[x], c i dcomplex[x], c i int64_t[x], c io dcomplex[xx], c io dcomplex[xx], c io dcomplex[xx])';
[C1, C2, C3] = AxiStokes3D_mex(mex_id_, nt, tx, nq, sx, M, C1, C2, C3, 1, nt, 1, nq, 1, nt, c3m, nt, c3m, nt, c3m);
C1=reshape(C1,nt,3*nq,Mp1); C2=reshape(C2,nt,3*nq,Mp1); C3=reshape(C3,nt,3*nq,Mp1);
end

% DLP (stresslet) pressure n-mode split coefs (scalar target): C1(log),C2(smooth),C3(Cauchy),C4(deriv-Cauchy); rows nt, no mu, uses source normal snx
