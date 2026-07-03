function A = axps_closeasm_mex(iphys, ilayer, nc, nt, tx, t3dx, p, np, nang, s3dx, s3dnx, s3dw, mu, ntcx, tcxi, idxall, S_ij, iself, A)
iphys=double(iphys); ilayer=double(ilayer); nc=double(nc); nt=double(nt); p=double(p); np=double(np);
nang=double(nang); ntcx=double(ntcx); iself=double(iself); mu=double(mu);
tx=double(tx(:)); tcxi=double(tcxi(:)); idxall=double(idxall(:)); s3dw=double(s3dw(:));
nsp=p*np; npp1=np+1; nsf=nsp*nang; nt3=3*nt; ns3=3*nsf; nap=nang*p;
t3dx=double(reshape(t3dx,3,nt)); s3dx=double(reshape(s3dx,3,nsf)); s3dnx=double(reshape(s3dnx,3,nsf));
S_ij=double(reshape(S_ij,ntcx,nap));
if iself, nrA=nsf; else, nrA=nt; end
A=double(reshape(A,nrA,nsf));
mex_id_ = 'axps_closeasm_r64(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i dcomplex[x], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i double[x], c i double[x], c i double[x], c i int64_t[x], c i double[x], c i double[x], c i double[xx], c i int64_t[x], c io double[xx])';
[A] = AxiStokes3D_mex(mex_id_, iphys, ilayer, nc, nt, tx, t3dx, p, np, nang, s3dx, s3dnx, s3dw, mu, ntcx, tcxi, idxall, S_ij, iself, A, 1, 1, 1, 1, nt, nt3, 1, 1, 1, ns3, ns3, nsf, 1, 1, npp1, ntcx, ntcx, nap, 1, nrA, nsf);
end

% ==== physical-space SPARSE close-correction APPLY (matvec piece; kernel-free) ====
