function K = axa_kernel_dlp_mex(nt, ns, srcr, srcz, srcnr, srcnz, tgtr, tgtz, m, mu, K)
nt=double(nt); ns=double(ns); m=double(m); mu=double(mu);
srcr=double(srcr(:)); srcz=double(srcz(:)); srcnr=double(srcnr(:)); srcnz=double(srcnz(:));
tgtr=double(tgtr(:)); tgtz=double(tgtz(:));
n3t=3*nt; n3s=3*ns; K=zeros(n3t,n3s);
mex_id_ = 'axa_kernel_dlp_r64(c i int64_t[x], c i int64_t[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i int64_t[x], c i double[x], c io double[xx])';
[K] = AxiStokes3D_mex(mex_id_, nt, ns, srcr, srcz, srcnr, srcnz, tgtr, tgtz, m, mu, K, 1, 1, ns, ns, ns, ns, nt, nt, 1, 1, n3t, n3s);
end

% Naive modal SLPn kernel (port of axissymstok_kernel_slpn.m).  TARGET normal (tgtnr,tgtnz,tgtnth);
% returns the meridional block K and the swirl block Ksw (Ksw prop. tgtnth; pass tgtnth=ones for unit).
