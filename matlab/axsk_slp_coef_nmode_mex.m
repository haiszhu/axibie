function [C1, C2, C3] = axsk_slp_coef_nmode_mex(nt, tx, nq, sx, M, mu)
nt=double(nt); nq=double(nq); M=double(M); mu=double(mu); tx=double(tx(:)); sx=double(sx(:));
Mp1=M+1; nt3=3*nt; c3m=3*nq*Mp1;
C1=complex(zeros(nt3,c3m)); C2=complex(zeros(nt3,c3m)); C3=complex(zeros(nt3,c3m));
mex_id_ = 'axsk_slp_coef_nmode_r64(c i int64_t[x], c i dcomplex[x], c i int64_t[x], c i dcomplex[x], c i int64_t[x], c i double[x], c io dcomplex[xx], c io dcomplex[xx], c io dcomplex[xx])';
[C1, C2, C3] = AxiStokes3D_mex(mex_id_, nt, tx, nq, sx, M, mu, C1, C2, C3, 1, nt, 1, nq, 1, 1, nt3, c3m, nt3, c3m, nt3, c3m);
C1=reshape(C1,nt3,3*nq,Mp1); C2=reshape(C2,nt3,3*nq,Mp1); C3=reshape(C3,nt3,3*nq,Mp1);
end

% SLPn n-mode split coefs (target normal tnx): C1,C2,C3,C4
