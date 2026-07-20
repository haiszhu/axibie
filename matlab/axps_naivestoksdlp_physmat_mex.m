function [As3d, Ad3d] = axps_naivestoksdlp_physmat_mex(nt, tx, t3dx, p, np, nang, sx, snx, sws, pmodes, mu, ifself, nrA, As3d, Ad3d)
nt=double(nt); p=double(p); np=double(np); nang=double(nang); M=double(pmodes); mu=double(mu(1)); ifself=double(ifself(1)); nrA=double(nrA(1));
tx=double(tx(:)); t3dx=double(t3dx); sx=double(sx(:)); N=np*p; Nnang3=3*N*nang;
if isempty(snx), snx=complex(zeros(N,1)); else, snx=complex(double(snx(:))); end
if isempty(sws), sws=zeros(N,1); else, sws=double(sws(:)); end
As3d=zeros(nrA,Nnang3); Ad3d=zeros(nrA,Nnang3);
mex_id_ = 'axps_naivestoksdlp_physmat_r64(c i int64_t[x], c i dcomplex[x], c i double[xx], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i double[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c io double[xx], c io double[xx])';
[As3d, Ad3d] = AxiStokes3D_mex(mex_id_, nt, tx, t3dx, p, np, nang, sx, snx, sws, M, mu, ifself, nrA, As3d, Ad3d, 1, nt, 3, nt, 1, 1, 1, N, N, N, 1, 1, 1, 1, nrA, Nnang3, nrA, Nnang3);
end

% whole-particle naive (smooth-quadrature) Stokes S'+D' (TARGET traction, sigma.n) physical operators (LAB xyz innermost)
% t3dnx = lab-frame 3D target normals (3 x nt), REQUIRED
% ifself=1: [3*N*nang x 3*N*nang] SELF (theta=0 + rotz tile); ifself=0: [3*nt x 3*N*nang] field
