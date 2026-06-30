function Mat = axt_lagrange_interp_mex(ns, xs, nt, xt)
ns=double(ns); nt=double(nt); xs=double(xs(:)); xt=double(xt(:));
Mat=zeros(nt,ns);
mex_id_ = 'axt_lagrange_interp_r64(c i int64_t[x], c i double[x], c i int64_t[x], c i double[x], c io double[xx])';
[Mat] = AxiStokes3D_mex(mex_id_, ns, xs, nt, xt, Mat, 1, ns, 1, nt, nt, ns);
end

% Real matrix × complex vector in Fortran arithmetic: y(m) = A(m,n) * x(n).
% Matches the block builders' arithmetic (-ffp-contract=off, no FMA), unlike MATLAB BLAS.
