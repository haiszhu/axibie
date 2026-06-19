function [K,Ksw] = axa_kernel_slpn_mex(nt, ns, srcr, srcz, tgtr, tgtz, tgtnr, tgtnz, tgtnth, m, mu, K, Ksw)
nt=double(nt); ns=double(ns); m=double(m); mu=double(mu);
srcr=double(srcr(:)); srcz=double(srcz(:)); tgtr=double(tgtr(:)); tgtz=double(tgtz(:));
tgtnr=double(tgtnr(:)); tgtnz=double(tgtnz(:)); tgtnth=double(tgtnth(:));
n3t=3*nt; n3s=3*ns; K=zeros(n3t,n3s); Ksw=zeros(n3t,n3s);
mex_id_ = 'axa_kernel_slpn_r64(c i int64_t[x], c i int64_t[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i double[x], c i int64_t[x], c i double[x], c io double[xx], c io double[xx])';
[K, Ksw] = AxiStokes3D_mex(mex_id_, nt, ns, srcr, srcz, tgtr, tgtz, tgtnr, tgtnz, tgtnth, m, mu, K, Ksw, 1, 1, ns, ns, nt, nt, nt, nt, nt, 1, 1, n3t, n3s, n3t, n3s);
end

% ------------------------------------------------------------
% SLP per-panel close-eval drivers (axissymstok_kernelsplit_mod; carrierVn/carrierFn/slp9coef
% two-carrier pipeline).  Complex inputs split real/imag.  Host supplies the upsampled panel:
% source nodes y=Z(tt) (nq=2p), source unit normals nv, real GL weights ws, complex speed
% weights wxp, panel endpoints za=Z(a0)/zb=Z(b0), iside (1=ext,0=int).
% --------------------------------------------------------------------------
% n-mode: A is the (nt,nq,9) complex block; returned flat as (nt, nq*9), reshape(A,nt,nq,9).
