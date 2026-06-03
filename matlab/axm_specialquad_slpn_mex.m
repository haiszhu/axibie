function [W,Wsw] = axm_specialquad_slpn_mex(nt, p, M, mu, hm, iside, zr, zi, nxr, nxi, yr, yi, ypr, ypi, zar, zai, zbr, zbi, gw, Pmat, W, Wsw)
nt=double(nt); p=double(p); M=double(M); mu=double(mu); hm=double(hm); iside=double(iside);
zr=double(zr(:)); zi=double(zi(:)); nxr=double(nxr(:)); nxi=double(nxi(:));
yr=double(yr(:)); yi=double(yi(:)); ypr=double(ypr(:)); ypi=double(ypi(:));
zar=double(zar); zai=double(zai); zbr=double(zbr); zbi=double(zbi); gw=double(gw(:)); Pmat=double(Pmat);
q=2*p; Mp1=M+1; NN=9*nt*p;
W=zeros(Mp1,NN); Wsw=zeros(Mp1,NN);
mex_id_ = 'axm_specialquad_slpn_r64(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i double[x], c i int64_t[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[xx], c io double[xx], c io double[xx])';
[W, Wsw] = AxiStokes3D_mex(mex_id_, nt, p, M, mu, hm, iside, zr, zi, nxr, nxi, yr, yi, ypr, ypi, zar, zai, zbr, zbi, gw, Pmat, W, Wsw, 1, 1, 1, 1, 1, 1, nt, nt, nt, nt, q, q, q, q, 1, 1, 1, 1, q, q, p, Mp1, NN, Mp1, NN);
end

% SLP close-eval operator block (port of StokAxiModalSpecialquad).  Same I/O shape as the DLP.
