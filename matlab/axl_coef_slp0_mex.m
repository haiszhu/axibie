function [C1,C2] = axl_coef_slp0_mex(nt, nq, zr, zi, yr, yi, C1, C2)
nt=double(nt); nq=double(nq);
zr=double(zr(:)); zi=double(zi(:)); yr=double(yr(:)); yi=double(yi(:));
C1=zeros(nt,nq); C2=zeros(nt,nq);
mex_id_ = 'axl_coef_slp0_r64(c i int64_t[x], c i int64_t[x], c i double[x], c i double[x], c i double[x], c i double[x], c io double[xx], c io double[xx])';
[C1, C2] = AxiStokes3D_mex(mex_id_, nt, nq, zr, zi, yr, yi, C1, C2, 1, 1, nt, nt, nq, nq, nt, nq, nt, nq);
end

% SLPn  S' : TARGET normal nxr,nxi.  C1 log, C2 smooth, C3 target-Cauchy density.
