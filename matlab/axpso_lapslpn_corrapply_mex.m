function u = axpso_lapslpn_corrapply_mex(K, x, u)
K=double(K); n=numel(x);
x=double(x(:)); u=double(u(:));
mex_id_ = 'axpso_lapslpn_corrapply_r64(c i int64_t[x], c i int64_t[x], c i double[x], c io double[x])';
[u] = AxiStokes3D_mex(mex_id_, K, n, x, u, 1, 1, n, n);
end

% unified handle-based near-correction setup: (ikernel,ilayer,params,iinter) x geometry -> handle
% ikernel 1 lap|2 stok|3 helm|4 maxw ; ilayer 1 SLP|2 SLPn|3 DLP|4 DLPn ; iinter 1 self|2 cross
% params = scalar kernel parameter (complex): wavenumber (helm/maxw); stokes reads real(params)
% as viscosity mu; laplace unused (0).
