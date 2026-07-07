function [Sk, ik, tcxik] = axpso_lapslpn_getself_mex(k, ntcxk, p, np, nang, Sk, ik, tcxik)
k=double(k); ntcxk=double(ntcxk); p=double(p); np=double(np); nang=double(nang);
npp1=np+1; r1=ntcxk; c1=nang*p;
Sk=double(reshape(Sk,r1,c1)); ik=double(ik(:)); tcxik=double(tcxik(:));
mex_id_ = 'axpso_lapslpn_getself_r64(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c io double[xx], c io double[x], c io double[x])';
[Sk, ik, tcxik] = AxiStokes3D_mex(mex_id_, k, ntcxk, p, np, nang, Sk, ik, tcxik, 1, 1, 1, 1, 1, r1, c1, ntcxk, npp1);
end

