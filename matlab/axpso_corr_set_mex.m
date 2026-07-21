function axpso_corr_set_mex(handle, vals, iks)
handle=double(handle); vals=double(vals(:)); iks=double(iks(:));
nvals=numel(vals); niks=numel(iks);
mex_id_ = 'axpso_corr_set_r64(c i double[x], c i int64_t[x], c i double[x], c i int64_t[x], c i double[x])';
AxiStokes3D_mex(mex_id_, handle, nvals, vals, niks, iks, 1, 1, nvals, 1, niks);
end

% assemble the stored correction into the caller's DENSE operator A (ACCUMULATES; seed A with the
% dense naive part) -- the explicit-matrix twin of axpso_corr_apply_mex.  Self handles: circulant +
% R(phi_a) conjugation + R_k sandwich; cross/eval handles: global canon target rows (rectangular ok).
