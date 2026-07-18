function [As3d, Ad3d] = axps_closelapsdlp_panel_mex(nt, tx, t3dx, tn3dx, p, nang, sx, snx, sws, swxp, sxlo, sxhi, tpan, gate, s3dx, s3dnx, s3dw, pmodes, iside, iclosed, ntcx, tcxi, idxall, iform, As3d, Ad3d)
nt=double(nt); p=double(p); nang=double(nang); M=double(pmodes); iside=double(iside); nangp=nang*p;
tx=double(tx(:)); t3dx=double(t3dx); sx=double(sx(:)); sxlo=double(sxlo(1)); sxhi=double(sxhi(1));
if isempty(tn3dx), tn3dx=zeros(3,nt); else, tn3dx=double(tn3dx); end
if isempty(snx), snx=complex(zeros(p,1)); else, snx=complex(double(snx(:))); end
if isempty(sws), sws=zeros(p,1); else, sws=double(sws(:)); end
if isempty(swxp), swxp=complex(zeros(p,1)); else, swxp=complex(double(swxp(:))); end
if isempty(tpan), tpan=zeros(2,1); else, tpan=double(tpan(:)); end
if isempty(gate), gate=0.0; else, gate=double(gate(1)); end
if isempty(s3dx), s3dx=zeros(3,1); else, s3dx=double(s3dx(:)); end
if isempty(s3dnx), s3dnx=zeros(3,1); else, s3dnx=double(s3dnx(:)); end
if isempty(s3dw), s3dw=zeros(1,1); else, s3dw=double(s3dw(1)); end
if isempty(iclosed), iclosed=0; else, iclosed=double(iclosed(1)); end
if isempty(ntcx), ntcx=0; else, ntcx=double(ntcx(1)); end
if isempty(tcxi), tcxi=zeros(1,1); else, tcxi=double(tcxi(:)); end
if isempty(idxall), idxall=zeros(1,1); else, idxall=double(idxall(:)); end
if isempty(iform), iform=0; else, iform=double(iform(1)); end
As3d=zeros(nt,nangp); Ad3d=zeros(nt,nangp);
mex_id_ = 'axps_closelapsdlp_panel_r64(c i int64_t[x], c i dcomplex[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i double[x], c i dcomplex[x], c i dcomplex[x], c i dcomplex[x], c i double[x], c i double[x], c i double[xx], c i double[xx], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i double[x], c i int64_t[x], c io double[xx], c io double[xx])';
[As3d, Ad3d] = AxiStokes3D_mex(mex_id_, nt, tx, t3dx, tn3dx, p, nang, sx, snx, sws, swxp, sxlo, sxhi, tpan, gate, s3dx, s3dnx, s3dw, M, iside, iclosed, ntcx, tcxi, idxall, iform, As3d, Ad3d, 1, nt, 3, nt, 3, nt, 1, 1, p, p, p, p, 1, 1, 2, 1, 3, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, nt, nangp, nt, nangp);
end

% single-panel Stokes S+D LP close-eval physical operators As3d, Ad3d [3*nt x 3*nang*p]
