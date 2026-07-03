function [Sk, ik, tcxik] = axpso_stokslpn_getself_mex(k, ntcxk, p, np, nang, Sk, ik, tcxik)
k=double(k); ntcxk=double(ntcxk); p=double(p); np=double(np); nang=double(nang);
npp1=np+1; r3=3*ntcxk; c3=3*nang*p;
Sk=double(reshape(Sk,r3,c3)); ik=double(ik(:)); tcxik=double(tcxik(:));
mex_id_ = 'axpso_stokslpn_getself_r64(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c io double[xx], c io double[x], c io double[x])';
[Sk, ik, tcxik] = AxiStokes3D_mex(mex_id_, k, ntcxk, p, np, nang, Sk, ik, tcxik, 1, 1, 1, 1, 1, r3, c3, ntcxk, npp1);
end

