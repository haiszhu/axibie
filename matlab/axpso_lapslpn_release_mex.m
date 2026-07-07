function ok = axpso_lapslpn_release_mex(ok)
ok=double(ok);
mex_id_ = 'axpso_lapslpn_release_r64(c io double[x])';
[ok] = AxiStokes3D_mex(mex_id_, ok, 1);
end

