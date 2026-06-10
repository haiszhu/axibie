function [K,Ksw] = axa_kernel_slpn_mex(nt, ns, srcr, srcz, tgtr, tgtz, tgtnr, tgtnz, tgtnth, m, mu, K, Ksw)
nt=double(nt); ns=double(ns); m=double(m); mu=double(mu);
srcr=double(srcr(:)); srcz=double(srcz(:)); tgtr=double(tgtr(:)); tgtz=double(tgtz(:));
tgtnr=double(tgtnr(:)); tgtnz=double(tgtnz(:)); tgtnth=double(tgtnth(:));
n3t=3*nt; n3s=3*ns; K=zeros(n3t,n3s); Ksw=zeros(n3t,n3s);
mex_id_ = 'axa_kernel_slpn_r64(c i int64_t[x], c i int64_t[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i int64_t[x], c i double[x], c io double[xx], c io double[xx])';
[K, Ksw] = AxiStokes3D_mex(mex_id_, nt, ns, srcr, srcz, tgtr, tgtz, tgtnr, tgtnz, tgtnth, m, mu, K, Ksw, 1, 1, ns, ns, nt, nt, nt, nt, nt, 1, 1, n3t, n3s, n3t, n3s);
end

% ============================================================
% Laplace 0th-mode kernel-split coefficients (axissymlap_kernelsplit_mod).
% All densities REAL, nt x nq.  The normals nu=t.nx, nu'=s.nx live in the special-quad
% matrices (Sspecialquad/Dspecialquad), NOT in these coefficients.  Naming follows the
% Stokes axa_coef_* wrappers with the 0th-mode "0" marker.
% ------------------------------------------------------------
% SLP  S : C1 log, C2 smooth.
