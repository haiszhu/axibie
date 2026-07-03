function [Sk, idxk, tcxik, canonk] = axpso_stokslpn_getcross_mex(k, ntcxk, nuk, p, np, nang, Sk, idxk, tcxik, canonk)
k=double(k); ntcxk=double(ntcxk); nuk=double(nuk); p=double(p); np=double(np); nang=double(nang);
npp1=np+1; r3=3*ntcxk; c3=3*nang*p;
Sk=double(reshape(Sk,r3,c3)); idxk=double(idxk(:)); tcxik=double(tcxik(:)); canonk=double(canonk(:));
mex_id_ = 'axpso_stokslpn_getcross_r64(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c io double[xx], c io double[x], c io double[x], c io double[x])';
[Sk, idxk, tcxik, canonk] = AxiStokes3D_mex(mex_id_, k, ntcxk, nuk, p, np, nang, Sk, idxk, tcxik, canonk, 1, 1, 1, 1, 1, 1, r3, c3, ntcxk, npp1, nuk);
end

