function [handle, nbytes] = axpso_corr_setup_mex(ikernel, ilayer, ipanel, dbg, params, iinter, K, p, np, pmodes, iside, iclosed, gate, geomoff, tpanoff, nsx, ntpan, sx, snx, sws, swxp, tpan, rball, ntarg, targoff, targ, targnx, Rm, Cc, handle)
ikernel=double(ikernel); ilayer=double(ilayer); ipanel=double(ipanel(1)); dbg=double(dbg(1)); iinter=double(iinter);
K=double(K); iside=double(iside); iclosed=double(iclosed);
params=complex(double(params)); gate=double(gate); rball=double(rball); ntarg=double(ntarg);
nsx=double(nsx); ntpan=double(ntpan); kp1=K+1; k9=9*K; k3=3*K;
p=double(p(:)); np=double(np(:)); pmodes=double(pmodes(:));
geomoff=double(geomoff(:)); tpanoff=double(tpanoff(:)); targoff=double(targoff(:));
sx=reshape(double(sx),nsx,1); snx=reshape(double(snx),nsx,1); swxp=reshape(double(swxp),nsx,1);
sws=reshape(double(sws),nsx,1); tpan=reshape(double(tpan),ntpan,1);
targ=double(reshape(targ,3,ntarg)); targnx=double(reshape(targnx,3,ntarg));
Rm=double(reshape(Rm,9,K)); Cc=double(reshape(Cc,3,K)); handle=double(handle); nbytes=0;
mex_id_ = 'axpso_corr_setup_r64(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i double[x], c i dcomplex[x], c i double[x], c i double[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[x], c i double[x], c io double[x], c io double[x])';
[handle, nbytes] = AxiStokes3D_mex(mex_id_, ikernel, ilayer, ipanel, dbg, params, iinter, K, p, np, pmodes, iside, iclosed, gate, geomoff, tpanoff, nsx, ntpan, sx, snx, sws, swxp, tpan, rball, ntarg, targoff, targ, targnx, Rm, Cc, handle, nbytes, 1, 1, 1, 1, 1, 1, 1, K, K, K, 1, 1, 1, kp1, kp1, 1, 1, nsx, nsx, nsx, nsx, ntpan, 1, 1, kp1, 3, ntarg, 3, ntarg, k9, k3, 1, 1);
end

% unified handle-based FULL close-eval operator setup -- SAME signature/workflow as corr_setup,
% but the fill stores the full close-eval value (for explicit assembly) instead of the correction.
% CROSS/EVAL targets are rotated into each source frame INSIDE the module (loc = R_k' (targ-C_k)).
