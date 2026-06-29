function u = axt_sto3dslp_eval_mex(tx, sx, sw, f)
nt=size(tx,2); ns=size(sx,2);
tx=double(tx); sx=double(sx); sw=double(sw(:)); f=double(f);
u=zeros(3,nt);
mex_id_ = 'axt_sto3dslp_eval_r64(c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i double[xx], c io double[xx])';
[u] = AxiStokes3D_mex(mex_id_, nt, tx, ns, sx, sw, f, u, 1, 3, nt, 1, 3, ns, ns, ns, 3, 3, nt);
end

% matrix-free 3D Stokes DLP (stresslet) velocity: u(3,nt) at targets tx(3,nt) from
% sources sx(3,ns) with normals snx(3,ns), weights sw(ns), force density f(ns,3)=[fx fy fz].
