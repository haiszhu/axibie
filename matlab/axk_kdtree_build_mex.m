function axk_kdtree_build_mex(backend, pts, leafsize)
backend=double(backend); leafsize=double(leafsize); pts=double(pts);
np=size(pts,1); ptsT=pts.';   % 3-by-np col-major == interleaved xyz for Fortran pts(3,np)
mex_id_ = 'axk_kdtree_build_r64(c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x])';
AxiStokes3D_mex(mex_id_, backend, np, ptsT, leafsize, 1, 1, 3, np, 1);
end

% all points within `radius` of qpt(3); returns 1-based indices sliced to the count found
