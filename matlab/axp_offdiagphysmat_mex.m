function [B, F, near] = axp_offdiagphysmat_mex(iphys, ilayer, nc, nt, tx, tphi, Ct, Rt, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, Cs, Rs, M, iside, iclosed, mu, iforce, tnx, tnphi)
iphys=double(iphys); ilayer=double(ilayer); nc=double(nc); nt=double(nt); p=double(p); np=double(np); M=double(M); iside=double(iside); iclosed=double(iclosed); mu=double(mu); iforce=double(iforce);
if nargin < 25 || isempty(tnx),   tnx   = zeros(nt,1); end             % target normal (meridian rho+i z, target frame): slpn/dlpn only
if nargin < 26 || isempty(tnphi), tnphi = zeros(nt,1); end             % target azimuthal normal component (target frame)
tx=double(tx(:)); tphi=double(tphi(:)); tnx=complex(double(tnx(:))); tnphi=double(tnphi(:));
sx=double(sx(:)); snx=double(snx(:)); swxp=double(swxp(:)); sxlo=double(sxlo(:)); sxhi=double(sxhi(:));
sws=double(sws(:)); tpan=double(tpan(:)); Ct=double(Ct(:)); Cs=double(Cs(:)); Rt=double(Rt(:)); Rs=double(Rs(:));
nsp=p*np; npp1=np+1; nphi=2*M+1; ncsp=nc*nsp; Nsrc=ncsp*nphi; Np1=Nsrc+1; nnzF=nphi*nphi*ncsp; cb=ncsp*nphi; nr=nc*nt;
B=complex(zeros(nr,cb)); iF=zeros(Np1,1); jF=zeros(nnzF,1); xF=complex(zeros(nnzF,1)); near=zeros(nt,1);
mex_id_ = 'axp_offdiagphysmat_r64(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i dcomplex[x], c i double[x], c i dcomplex[x], c i double[x], c i double[x], c i double[x], c i int64_t[x], c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i double[x], c i dcomplex[x], c i double[x], c i dcomplex[x], c i dcomplex[x], c i double[x], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io dcomplex[xx], c io double[x], c io double[x], c io dcomplex[x], c io double[x])';
[B, iF, jF, xF, near] = AxiStokes3D_mex(mex_id_, iphys, ilayer, nc, nt, tx, tphi, tnx, tnphi, Ct, Rt, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, Cs, Rs, M, iside, iclosed, mu, iforce, B, iF, jF, xF, near, 1, 1, 1, 1, nt, nt, nt, nt, 3, 9, 1, 1, nsp, nsp, nsp, nsp, npp1, np, np, 3, 9, 1, 1, 1, 1, 1, nr, cb, Np1, nnzF, nnzF, nt);
F = sparse(repelem((1:Nsrc)', diff(round(iF))), jF, xF, Nsrc, Nsrc);
near = near > 0.5;
end

% ============================================================
% Helpers for the pure-MATLAB n-mode block builders (utils/axls_*_blockmat_nmode.m):
% the SAME compiled gauss / lagrange-interp / carrier / sdspecialquad / split-coef routines.
% ============================================================

% Gauss-Legendre nodes/weights + spectral diff matrix on [-1,1]
