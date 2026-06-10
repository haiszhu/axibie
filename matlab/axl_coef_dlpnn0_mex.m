function [C1,C2,C3a,C3b,C4a,C4b,C5] = axl_coef_dlpnn0_mex(nt, nq, zr, zi, yr, yi, nxr, nxi, nvr, nvi, C1, C2, C3a, C3b, C4a, C4b, C5)
nt=double(nt); nq=double(nq);
zr=double(zr(:)); zi=double(zi(:)); yr=double(yr(:)); yi=double(yi(:));
nxr=double(nxr(:)); nxi=double(nxi(:)); nvr=double(nvr(:)); nvi=double(nvi(:));
C1=zeros(nt,nq); C2=zeros(nt,nq); C3a=zeros(nt,nq); C3b=zeros(nt,nq); C4a=zeros(nt,nq); C4b=zeros(nt,nq); C5=zeros(nt,nq);
mex_id_ = 'axl_coef_dlpnn0_r64(c i int64_t[x], c i int64_t[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c io double[xx], c io double[xx], c io double[xx], c io double[xx], c io double[xx], c io double[xx], c io double[xx])';
[C1, C2, C3a, C3b, C4a, C4b, C5] = AxiStokes3D_mex(mex_id_, nt, nq, zr, zi, yr, yi, nxr, nxi, nvr, nvi, C1, C2, C3a, C3b, C4a, C4b, C5, 1, 1, nt, nt, nq, nq, nt, nt, nq, nq, nt, nq, nt, nq, nt, nq, nt, nq, nt, nq, nt, nq, nt, nq);
end

% ============================================================
% Laplace 0th-mode close-eval OPERATORS (axissymlap_specialquad_mod).
% Caller upsamples the source to q=2p (y=Z(tt), yp=Zp(tt)) and supplies the q->p Legendre
% projection Pmat (q x p), GL weights gw (q), endpoints za=Z(a0)/zb=Z(b0), hm=(b0-a0)/2,
% iside (1=ext,0=int).  W is the nt x p close-eval block.  Complex args split real/imag.
% ------------------------------------------------------------
% SLP  S
