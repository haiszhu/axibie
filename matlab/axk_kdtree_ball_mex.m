function [idx, nout] = axk_kdtree_ball_mex(qpt, radius, nmax)
qpt=double(qpt(:)); radius=double(radius); nmax=double(nmax);
idx=zeros(nmax,1); nout=0;
mex_id_ = 'axk_kdtree_ball_r64(c i double[x], c i double[x], c i int64_t[x], c io int64_t[x], c io int64_t[x])';
[idx, nout] = AxiStokes3D_mex(mex_id_, qpt, radius, nmax, idx, nout, 3, 1, 1, nmax, 1);
idx=idx(1:min(nout,nmax))+1;   % 0-based -> 1-based
end

% release the held tree
