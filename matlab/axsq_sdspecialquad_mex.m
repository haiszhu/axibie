function [As, Ad, A1, A2, A3, A4] = axsq_sdspecialquad_mex(nt, zt, p, zsrc, nzsrc, wzp, za, zb, iside)
nt=double(nt); p=double(p); iside=double(iside);
zt=double(zt(:)); zsrc=double(zsrc(:)); nzsrc=double(nzsrc(:)); wzp=double(wzp(:));
za=double(za(1)); zb=double(zb(1));
As=zeros(p,nt); A1=zeros(p,nt); A2=zeros(p,nt); A3=zeros(p,nt); A4=zeros(p,nt); Ad=complex(zeros(p,nt));
mex_id_ = 'axsq_sdspecialquad_r64(c i int64_t[x], c i dcomplex[x], c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i dcomplex[x], c i dcomplex[x], c i dcomplex[x], c i int64_t[x], c io double[xx], c io dcomplex[xx], c io double[xx], c io double[xx], c io double[xx], c io double[xx])';
[As, Ad, A1, A2, A3, A4] = AxiStokes3D_mex(mex_id_, nt, zt, p, zsrc, nzsrc, wzp, za, zb, iside, As, Ad, A1, A2, A3, A4, 1, nt, 1, p, p, p, 1, 1, 1, p, nt, p, nt, p, nt, p, nt, p, nt, p, nt);
end

% single-panel Laplace S+D LP close-eval physical operators As3d, Ad3d [nt x nang*p]
