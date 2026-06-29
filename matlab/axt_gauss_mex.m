function [tgl, wgl, Dgl] = axt_gauss_mex(n)
n=double(n);
tgl=zeros(n,1); wgl=zeros(n,1); Dgl=zeros(n,n);
mex_id_ = 'axt_gauss_r64(c i int64_t[x], c io double[x], c io double[x], c io double[xx])';
[tgl, wgl, Dgl] = AxiStokes3D_mex(mex_id_, n, tgl, wgl, Dgl, 1, n, n, n, n);
end

% Lagrange interpolation matrix: values on nodes xs (ns) -> targets xt (nt)
