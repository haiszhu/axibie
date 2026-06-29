function [C1, C2, C3] = axlk_dlp_coef_nmode_mex(nt, tx, nq, sx, snx, M)
nt=double(nt); nq=double(nq); M=double(M); tx=double(tx(:)); sx=double(sx(:)); snx=double(snx(:)); Mp1=M+1; cm=nq*Mp1;
C1=zeros(nt,cm); C2=zeros(nt,cm); C3=zeros(nt,cm);
mex_id_ = 'axlk_dlp_coef_nmode_r64(c i int64_t[x], c i dcomplex[x], c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i int64_t[x], c io double[xx], c io double[xx], c io double[xx])';
[C1, C2, C3] = AxiStokes3D_mex(mex_id_, nt, tx, nq, sx, snx, M, C1, C2, C3, 1, nt, 1, nq, nq, 1, nt, cm, nt, cm, nt, cm);
C1=reshape(C1,nt,nq,Mp1); C2=reshape(C2,nt,nq,Mp1); C3=reshape(C3,nt,nq,Mp1);
end

% DLPn n-mode split coefs (both normals): C1,C2,C3a(target Cauchy),C3b(source Cauchy),C4(deriv-Cauchy)
