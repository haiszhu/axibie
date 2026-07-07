function u = axpso_stokslpn_corrapply_mex(K, Rm, x, u)
K=double(K); n=numel(x); k9=9*K;
Rm=double(reshape(Rm,9,K)); x=double(x(:)); u=double(u(:));
mex_id_ = 'axpso_stokslpn_corrapply_r64(c i int64_t[x], c i int64_t[x], c i double[x], c i double[x], c io double[x])';
[u] = AxiStokes3D_mex(mex_id_, K, n, Rm, x, u, 1, 1, k9, n, n);
end

