function A = axpso_corr2dense_get_mex(handle, nra, nca, A)
handle=double(handle); nra=double(nra); nca=double(nca);
A=double(A);
mex_id_ = 'axpso_corr2dense_get_r64(c i double[x], c i int64_t[x], c i int64_t[x], c io double[xx])';
[A] = AxiStokes3D_mex(mex_id_, handle, nra, nca, A, 1, 1, 1, nra, nca);
end

% Fortran port of the driver's close2corr: turn a FULL close-eval handle (axpso_close_setup_mex) into
% a CORRECTION handle in place, subtracting the naive Sto3d kernel in the module's own compact layout.
% Topology (idx/canon/tcxi/Rm, ilayer/mu/nc/interaction, p/np/nang) is read from the handle; only the
% geometry (source sx/snx/sws, targets targ/targnx, centres Cc) is passed.  Then corr_apply scatters it.
