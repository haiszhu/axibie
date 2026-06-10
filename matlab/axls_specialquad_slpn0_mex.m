function W = axls_specialquad_slpn0_mex(nt, p, hm, iside, zr, zi, nxr, nxi, yr, yi, ypr, ypi, zar, zai, zbr, zbi, gw, Pmat, W)
nt=double(nt); p=double(p); hm=double(hm); iside=double(iside);
zr=double(zr(:)); zi=double(zi(:)); nxr=double(nxr(:)); nxi=double(nxi(:));
yr=double(yr(:)); yi=double(yi(:)); ypr=double(ypr(:)); ypi=double(ypi(:));
zar=double(zar); zai=double(zai); zbr=double(zbr); zbi=double(zbi); gw=double(gw(:)); Pmat=double(Pmat);
q=2*p;
W=zeros(nt,p);
mex_id_ = 'axls_specialquad_slpn0_r64(c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[xx], c io double[xx])';
[W] = AxiStokes3D_mex(mex_id_, nt, p, hm, iside, zr, zi, nxr, nxi, yr, yi, ypr, ypi, zar, zai, zbr, zbi, gw, Pmat, W, 1, 1, 1, 1, nt, nt, nt, nt, q, q, q, q, 1, 1, 1, 1, q, q, p, nt, p);
end

% DLP  D : source normal derived from yp (no extra input).
