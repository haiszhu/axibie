function u = axpso_corr_apply_mex(handle, nx, x, nu, u)
handle=double(handle); nx=double(nx); nu=double(nu);
x=double(x(:)); u=double(u(:));
mex_id_ = 'axpso_corr_apply_r64(c i double[x], c i int64_t[x], c i double[x], c i int64_t[x], c io double[x])';
[u] = AxiStokes3D_mex(mex_id_, handle, nx, x, nu, u, 1, 1, nx, 1, nu);
end
