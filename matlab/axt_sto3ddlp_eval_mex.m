function u = axt_sto3ddlp_eval_mex(tx, sx, snx, sw, f)
nt=size(tx,2); ns=size(sx,2);
tx=double(tx); sx=double(sx); snx=double(snx); sw=double(sw(:)); f=double(f);
u=zeros(3,nt);
mex_id_ = 'axt_sto3ddlp_eval_r64(c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[xx], c i double[x], c i double[xx], c io double[xx])';
[u] = AxiStokes3D_mex(mex_id_, nt, tx, ns, sx, snx, sw, f, u, 1, 3, nt, 1, 3, ns, 3, ns, ns, ns, 3, 3, nt);
end

% All-modes Q-side carrier vk=Q_{n-1/2}(chi), ve, n=0..M  (scalar chi)
