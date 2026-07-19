function [As3d, Ad3d] = axps_naivelapsdlpn_physmat_mex(nt, tx, t3dx, t3dnx, p, np, nang, sx, snx, sws, pmodes, ifself, nrA, As3d, Ad3d)
nt=double(nt); p=double(p); np=double(np); nang=double(nang); M=double(pmodes); ifself=double(ifself(1)); nrA=double(nrA(1));
tx=double(tx(:)); t3dx=double(t3dx); t3dnx=double(t3dnx); sx=double(sx(:)); N=np*p; Nnang=N*nang;
if isempty(snx), snx=complex(zeros(N,1)); else, snx=complex(double(snx(:))); end
if isempty(sws), sws=zeros(N,1); else, sws=double(sws(:)); end
As3d=zeros(nrA,Nnang); Ad3d=zeros(nrA,Nnang);
mex_id_ = 'axps_naivelapsdlpn_physmat_r64(c i int64_t[x], c i dcomplex[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c io double[xx], c io double[xx])';
[As3d, Ad3d] = AxiStokes3D_mex(mex_id_, nt, tx, t3dx, t3dnx, p, np, nang, sx, snx, sws, M, ifself, nrA, As3d, Ad3d, 1, nt, 3, nt, 3, nt, 1, 1, 1, N, N, N, 1, 1, 1, nrA, Nnang, nrA, Nnang);
end

% SLP n-mode split coefs C1(log), C2(smooth)
