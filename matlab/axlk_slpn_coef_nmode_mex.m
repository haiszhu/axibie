function [C1, C2, C3] = axlk_slpn_coef_nmode_mex(nt, tx, tnx, nq, sx, M)
nt=double(nt); nq=double(nq); M=double(M); tx=double(tx(:)); tnx=double(tnx(:)); sx=double(sx(:)); Mp1=M+1; cm=nq*Mp1;
C1=zeros(nt,cm); C2=zeros(nt,cm); C3=zeros(nt,cm);
mex_id_ = 'axlk_slpn_coef_nmode_r64(c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i int64_t[x], c i dcomplex[x], c i int64_t[x], c io double[xx], c io double[xx], c io double[xx])';
[C1, C2, C3] = AxiStokes3D_mex(mex_id_, nt, tx, tnx, nq, sx, M, C1, C2, C3, 1, nt, nt, 1, nq, 1, nt, cm, nt, cm, nt, cm);
C1=reshape(C1,nt,nq,Mp1); C2=reshape(C2,nt,nq,Mp1); C3=reshape(C3,nt,nq,Mp1);
end

% DLP n-mode split coefs (source normal snx): C1,C2,C3(source Cauchy)
