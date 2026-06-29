function [C1, C2] = axlk_slp_coef_nmode_mex(nt, tx, nq, sx, M)
nt=double(nt); nq=double(nq); M=double(M); tx=double(tx(:)); sx=double(sx(:)); Mp1=M+1; cm=nq*Mp1;
C1=zeros(nt,cm); C2=zeros(nt,cm);
mex_id_ = 'axlk_slp_coef_nmode_r64(c i int64_t[x], c i dcomplex[x], c i int64_t[x], c i dcomplex[x], c i int64_t[x], c io double[xx], c io double[xx])';
[C1, C2] = AxiStokes3D_mex(mex_id_, nt, tx, nq, sx, M, C1, C2, 1, nt, 1, nq, 1, nt, cm, nt, cm);
C1=reshape(C1,nt,nq,Mp1); C2=reshape(C2,nt,nq,Mp1);
end

% SLPn n-mode split coefs (target normal tnx): C1,C2,C3(target Cauchy)
