function A = axp_physmat_setup_mex(ikernel, ilayer, params, iinter, K, p, np, pmodes, iside, iclosed, geomoff, tpanoff, nsx, ntpan, sx, snx, sws, swxp, tpan, ntarg, targoff, targ, targnx, Rm, Cc, nra, nca)
ikernel=double(ikernel); ilayer=double(ilayer); iinter=double(iinter);
K=double(K); iside=double(iside); iclosed=double(iclosed);
params=complex(double(params)); ntarg=double(ntarg);
nsx=double(nsx); ntpan=double(ntpan); kp1=K+1; k9=9*K; k3=3*K;
p=double(p(:)); np=double(np(:)); pmodes=double(pmodes(:));
geomoff=double(geomoff(:)); tpanoff=double(tpanoff(:)); targoff=double(targoff(:));
sx=reshape(double(sx),nsx,1); snx=reshape(double(snx),nsx,1); swxp=reshape(double(swxp),nsx,1);
sws=reshape(double(sws),nsx,1); tpan=reshape(double(tpan),ntpan,1);
targ=double(reshape(targ,3,ntarg)); targnx=double(reshape(targnx,3,ntarg));
Rm=double(reshape(Rm,9,K)); Cc=double(reshape(Cc,3,K));
nra=double(nra); nca=double(nca);
A=zeros(nra,nca);
mex_id_ = 'axp_physmat_setup_r64(c i int64_t[x], c i int64_t[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i double[x], c i dcomplex[x], c i double[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[x], c i double[x], c i int64_t[x], c i int64_t[x], c io double[xx])';
[A] = AxiStokes3D_mex(mex_id_, ikernel, ilayer, params, iinter, K, p, np, pmodes, iside, iclosed, geomoff, tpanoff, nsx, ntpan, sx, snx, sws, swxp, tpan, ntarg, targoff, targ, targnx, Rm, Cc, nra, nca, A, 1, 1, 1, 1, 1, K, K, K, 1, 1, kp1, kp1, 1, 1, nsx, nsx, nsx, nsx, ntpan, 1, kp1, 3, ntarg, 3, ntarg, k9, k3, 1, 1, nra, nca);
end

% ==== LEVEL-2 MODAL master (axissym_modemat_setup): K-particle ragged geometry -> per-particle
% all-mode block stacks A(nra, nca, max(pmodes)+1), complex for BOTH kernels.  iinter 1 self | 3 eval
% (2 cross fails loud in Fortran).  A is pre-zeroed HERE (the module writes only each particle's
% pmodes(k)+1 modes; higher modes of low-pmode particles stay 0).
