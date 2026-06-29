function [C1, C2, C3a, C3b, C4] = axlk_dlpn_coef_nmode_mex(nt, tx, tnx, nq, sx, snx, M)
nt=double(nt); nq=double(nq); M=double(M); tx=double(tx(:)); tnx=double(tnx(:)); sx=double(sx(:)); snx=double(snx(:)); Mp1=M+1; cm=nq*Mp1;
C1=zeros(nt,cm); C2=zeros(nt,cm); C3a=zeros(nt,cm); C3b=zeros(nt,cm); C4=zeros(nt,cm);
mex_id_ = 'axlk_dlpn_coef_nmode_r64(c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i int64_t[x], c io double[xx], c io double[xx], c io double[xx], c io double[xx], c io double[xx])';
[C1, C2, C3a, C3b, C4] = AxiStokes3D_mex(mex_id_, nt, tx, tnx, nq, sx, snx, M, C1, C2, C3a, C3b, C4, 1, nt, nt, 1, nq, nq, 1, nt, cm, nt, cm, nt, cm, nt, cm, nt, cm);
C1=reshape(C1,nt,nq,Mp1); C2=reshape(C2,nt,nq,Mp1); C3a=reshape(C3a,nt,nq,Mp1); C3b=reshape(C3b,nt,nq,Mp1); C4=reshape(C4,nt,nq,Mp1);
end
