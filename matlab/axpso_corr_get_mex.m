function [blk, ik, tcxik] = axpso_corr_get_mex(handle, k, ntcxk, nc, p, np, pmodes, blk, ik, tcxik)
handle=double(handle); k=double(k); ntcxk=double(ntcxk); nc=double(nc); p=double(p); np=double(np); pmodes=double(pmodes);
npp1=np+1; nb=nc*nc*ntcxk*(2*pmodes+1)*p;
blk=double(blk(:)); ik=double(ik(:)); tcxik=double(tcxik(:));
mex_id_ = 'axpso_corr_get_r64(c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c io double[x], c io double[x], c io double[x])';
[blk, ik, tcxik] = AxiStokes3D_mex(mex_id_, handle, k, ntcxk, nc, p, np, pmodes, blk, ik, tcxik, 1, 1, 1, 1, 1, 1, 1, nb, ntcxk, npp1);
end

% stage-3 SET (ONE call, whole record): deliver an externally computed record into the struct the
% handle points to.  vals = ALL value blocks flat (valoff layout); iks = pair map flat (off layout;
% self: azimuth-0 meridian index, cross: GLOBAL target id, canonicalized inside).  Lengths validated
% against the handle's topology.  Per-particle read stays axpso_corr_get_mex.
