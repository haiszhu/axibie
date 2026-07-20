function [Asn3d, Adn3d] = axps_naivestoksdlpn_physmat_mex(nt, tx, t3dx, t3dnx, p, np, nang, sx, snx, sws, pmodes, mu, ifself, nrA, Asn3d, Adn3d)
nt=double(nt); p=double(p); np=double(np); nang=double(nang); M=double(pmodes); mu=double(mu(1)); ifself=double(ifself(1)); nrA=double(nrA(1));
tx=double(tx(:)); t3dx=double(t3dx); t3dnx=double(t3dnx); sx=double(sx(:)); N=np*p; Nnang3=3*N*nang;
if isempty(snx), snx=complex(zeros(N,1)); else, snx=complex(double(snx(:))); end
if isempty(sws), sws=zeros(N,1); else, sws=double(sws(:)); end
Asn3d=zeros(nrA,Nnang3); Adn3d=zeros(nrA,Nnang3);
mex_id_ = 'axps_naivestoksdlpn_physmat_r64(c i int64_t[x], c i dcomplex[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i double[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c io double[xx], c io double[xx])';
[Asn3d, Adn3d] = AxiStokes3D_mex(mex_id_, nt, tx, t3dx, t3dnx, p, np, nang, sx, snx, sws, M, mu, ifself, nrA, Asn3d, Adn3d, 1, nt, 3, nt, 3, nt, 1, 1, 1, N, N, N, 1, 1, 1, 1, nrA, Nnang3, nrA, Nnang3);
end

% SLP n-mode split coefs C1(log), C2(smooth)
