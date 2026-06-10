function [C1,C2,C3] = axl_coef_slpn0_mex(nt, nq, zr, zi, yr, yi, nxr, nxi, C1, C2, C3)
nt=double(nt); nq=double(nq);
zr=double(zr(:)); zi=double(zi(:)); yr=double(yr(:)); yi=double(yi(:)); nxr=double(nxr(:)); nxi=double(nxi(:));
C1=zeros(nt,nq); C2=zeros(nt,nq); C3=zeros(nt,nq);
mex_id_ = 'axl_coef_slpn0_r64(c i int64_t[x], c i int64_t[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c io double[xx], c io double[xx], c io double[xx])';
[C1, C2, C3] = AxiStokes3D_mex(mex_id_, nt, nq, zr, zi, yr, yi, nxr, nxi, C1, C2, C3, 1, 1, nt, nt, nq, nq, nt, nt, nt, nq, nt, nq, nt, nq);
end

% DLP  D : SOURCE normal nvr,nvi.  C1 log, C2 smooth, C3 source-Cauchy density.
