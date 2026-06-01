function K = axa_kernel_slp_mex(nt, ns, srcr, srcz, tgtr, tgtz, m, mu, K)
nt=double(nt); ns=double(ns); m=double(m); mu=double(mu);
srcr=double(srcr(:)); srcz=double(srcz(:)); tgtr=double(tgtr(:)); tgtz=double(tgtz(:));
n3t=3*nt; n3s=3*ns; K=zeros(n3t,n3s);
mex_id_ = 'axa_kernel_slp_r64(c i int64_t[x], c i int64_t[x], c i double[x], c i double[x], c i double[x], c i double[x], c i int64_t[x], c i double[x], c io double[xx])';
[K] = AxiStokes3D_mex(mex_id_, nt, ns, srcr, srcz, tgtr, tgtz, m, mu, K, 1, 1, ns, ns, nt, nt, 1, 1, n3t, n3s);
end

% Naive modal DLP kernel (port of axissymstok_kernel_dlp.m).  Needs source normals srcnr,srcnz.
