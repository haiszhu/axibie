function y = axt_rmatcz_mex(m, n, A, x)
m=double(m); n=double(n); A=double(A); x=double(x(:));
y=complex(zeros(m,1));
mex_id_ = 'axt_rmatcz_r64(c i int64_t[x], c i int64_t[x], c i double[xx], c i dcomplex[x], c io dcomplex[x])';
[y] = AxiStokes3D_mex(mex_id_, m, n, A, x, y, 1, 1, m, n, n, m);
end

% Evaluate Lagrange interpolant yout(j)=sum_k Mat(j,k)*yc(k) in Fortran arithmetic.
% xs(ns) nodes, yc(ns) complex values at nodes, xt(nt) eval points → yout(nt) complex.
% matmul done in Fortran so arithmetic matches the block builders' internal ec2=matmul(IPe2,xc).
