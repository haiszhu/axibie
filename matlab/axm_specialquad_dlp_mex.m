function W = axm_specialquad_dlp_mex(nt, p, M, mu, hm, iside, zr, zi, yr, yi, ypr, ypi, zar, zai, zbr, zbi, gw, Pmat, W)
nt=double(nt); p=double(p); M=double(M); mu=double(mu); hm=double(hm); iside=double(iside);
zr=double(zr(:)); zi=double(zi(:)); yr=double(yr(:)); yi=double(yi(:)); ypr=double(ypr(:)); ypi=double(ypi(:));
zar=double(zar); zai=double(zai); zbr=double(zbr); zbi=double(zbi); gw=double(gw(:)); Pmat=double(Pmat);
q=2*p; Mp1=M+1; NN=9*nt*p;
W=zeros(Mp1,NN);
mex_id_ = 'axm_specialquad_dlp_r64(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i double[x], c i int64_t[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[xx], c io double[xx])';
[W] = AxiStokes3D_mex(mex_id_, nt, p, M, mu, hm, iside, zr, zi, yr, yi, ypr, ypi, zar, zai, zbr, zbi, gw, Pmat, W, 1, 1, 1, 1, 1, 1, nt, nt, q, q, q, q, 1, 1, 1, 1, q, q, p, Mp1, NN);
end

% SLPn close-eval operator block (port of StokSLPnAxiModalSpecialquad_v2).  TARGET normal nx
% (length nt) split nxr/nxi.  TWO-BLOCK: meridional W and swirl Wsw, each (M+1)x(9 nt p);
% reshape to (M+1,3nt,3p).  Caller scales Wsw by n_theta + applies the flipped parity mask.
