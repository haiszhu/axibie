function axpso_close2corr_set_mex(handle, K, geomoff, nsx, sx, snx, sws, targoff, ntarg, targ, targnx, Cc)
handle=double(handle); K=double(K); nsx=double(nsx); ntarg=double(ntarg); kp1=K+1; k3=3*K;
geomoff=double(geomoff(:)); targoff=double(targoff(:));
sx=reshape(double(sx),nsx,1); snx=reshape(double(snx),nsx,1); sws=reshape(double(sws),nsx,1);
targ=double(reshape(targ,3,ntarg)); targnx=double(reshape(targnx,3,ntarg)); Cc=double(reshape(Cc,3,K));
mex_id_ = 'axpso_close2corr_set_r64(c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i double[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[x])';
AxiStokes3D_mex(mex_id_, handle, K, geomoff, nsx, sx, snx, sws, targoff, ntarg, targ, targnx, Cc, 1, 1, kp1, 1, nsx, nsx, nsx, kp1, 1, 3, ntarg, 3, ntarg, k3);
end

% unified handle-based apply: source density x (nx) -> target buffer u (nu), ACCUMULATES (seed u
% with the naive FMM part).  ONE apply spans self / cross / arbitrary-target eval; for self and
% cross nu==nx, for arbitrary targets nu = target space.  Kernel/layer/interaction are in the handle.
