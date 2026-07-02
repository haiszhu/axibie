function [tcxi, ntcx] = axps_closesize_mex(nt, tx, t3dx, p, np, sx, sws, gate, tcxi, ntcx)
nt=double(nt); p=double(p); np=double(np); gate=double(gate); ntcx=double(ntcx);
tx=double(tx(:)); sx=double(sx(:)); sws=double(sws(:)); tcxi=double(tcxi(:));
nsp=p*np; npp1=np+1; nt3=3*nt;
t3dx=double(reshape(t3dx,3,nt));
mex_id_ = 'axps_closesize_r64(c i int64_t[x], c i dcomplex[x], c i double[x], c i int64_t[x], c i int64_t[x], c i dcomplex[x], c i double[x], c i double[x], c io double[x], c io double[x])';
[tcxi, ntcx] = AxiStokes3D_mex(mex_id_, nt, tx, t3dx, p, np, sx, sws, gate, tcxi, ntcx, 1, nt, nt3, 1, 1, nsp, nsp, 1, npp1, 1);
end

% ==== physical-space SPARSE close-correction (pass 3/3): ASSEMBLE naive + scatter ====
