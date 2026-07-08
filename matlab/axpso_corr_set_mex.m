function axpso_corr_set_mex(handle, k, ntcxk, nc, p, np, pmodes, blk)
handle=double(handle); k=double(k); ntcxk=double(ntcxk); nc=double(nc); p=double(p); np=double(np); pmodes=double(pmodes);
nb=nc*nc*ntcxk*(2*pmodes+1)*p;
blk=double(blk(:));
mex_id_ = 'axpso_corr_set_r64(c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x])';
AxiStokes3D_mex(mex_id_, handle, k, ntcxk, nc, p, np, pmodes, blk, 1, 1, 1, 1, 1, 1, 1, nb);
end

% assemble the stored correction into the caller's DENSE operator A (ACCUMULATES; seed A with the
% dense naive part) -- the explicit-matrix twin of axpso_corr_apply_mex.  Self handles: circulant +
% R(phi_a) conjugation + R_k sandwich; cross/eval handles: global canon target rows (rectangular ok).
