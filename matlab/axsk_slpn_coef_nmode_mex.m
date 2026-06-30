function [C1, C2, C3, C4] = axsk_slpn_coef_nmode_mex(nt, tx, tnx, nq, sx, M, mu)
nt=double(nt); nq=double(nq); M=double(M); mu=double(mu); tx=double(tx(:)); tnx=double(tnx(:)); sx=double(sx(:));
Mp1=M+1; nt3=3*nt; c3m=3*nq*Mp1;
C1=complex(zeros(nt3,c3m)); C2=complex(zeros(nt3,c3m)); C3=complex(zeros(nt3,c3m)); C4=complex(zeros(nt3,c3m));
mex_id_ = 'axsk_slpn_coef_nmode_r64(c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i int64_t[x], c i dcomplex[x], c i int64_t[x], c i double[x], c io dcomplex[xx], c io dcomplex[xx], c io dcomplex[xx], c io dcomplex[xx])';
[C1, C2, C3, C4] = AxiStokes3D_mex(mex_id_, nt, tx, tnx, nq, sx, M, mu, C1, C2, C3, C4, 1, nt, nt, 1, nq, 1, 1, nt3, c3m, nt3, c3m, nt3, c3m, nt3, c3m);
C1=reshape(C1,nt3,3*nq,Mp1); C2=reshape(C2,nt3,3*nq,Mp1); C3=reshape(C3,nt3,3*nq,Mp1); C4=reshape(C4,nt3,3*nq,Mp1);
end

% DLP n-mode split coefs (source normal snx): C1,C2,C3,C4
