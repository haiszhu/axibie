function [C1,C2,C3a,C3b,C4] = axl_coef_dlpn0_mex(nt, nq, zr, zi, yr, yi, nxr, nxi, nvr, nvi, C1, C2, C3a, C3b, C4)
nt=double(nt); nq=double(nq);
zr=double(zr(:)); zi=double(zi(:)); yr=double(yr(:)); yi=double(yi(:));
nxr=double(nxr(:)); nxi=double(nxi(:)); nvr=double(nvr(:)); nvi=double(nvi(:));
C1=zeros(nt,nq); C2=zeros(nt,nq); C3a=zeros(nt,nq); C3b=zeros(nt,nq); C4=zeros(nt,nq);
mex_id_ = 'axl_coef_dlpn0_r64(c i int64_t[x], c i int64_t[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c io double[xx], c io double[xx], c io double[xx], c io double[xx], c io double[xx])';
[C1, C2, C3a, C3b, C4] = AxiStokes3D_mex(mex_id_, nt, nq, zr, zi, yr, yi, nxr, nxi, nvr, nvi, C1, C2, C3a, C3b, C4, 1, 1, nt, nt, nq, nq, nt, nt, nq, nq, nt, nq, nt, nq, nt, nq, nt, nq, nt, nq);
end

% DLPnn  D'' : TARGET + SOURCE normals.  C1 log, C2 smooth, C3a target-Cauchy,
% C3b source-Cauchy, C4a (nu^2) deriv-Cauchy, C4b (nu nu') deriv-Cauchy, C5 cubic.
