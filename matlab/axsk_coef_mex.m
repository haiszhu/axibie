function [C1, C2, C3, C4, C5] = axsk_coef_mex(ilayer, nt, tx, tnx, nq, sx, snx, M, mu)
ilayer=double(ilayer); nt=double(nt); nq=double(nq); M=double(M); mu=double(mu); tx=double(tx(:)); tnx=double(tnx(:)); sx=double(sx(:)); snx=double(snx(:));
Mp1=M+1; if ilayer<=4, nrow=3*nt; else, nrow=nt; end; c3m=3*nq*Mp1;
C1=complex(zeros(nrow,c3m)); C2=complex(zeros(nrow,c3m)); C3=complex(zeros(nrow,c3m)); C4=complex(zeros(nrow,c3m)); C5=complex(zeros(nrow,c3m));
mex_id_ = 'axsk_coef_r64(c i int64_t[x], c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i int64_t[x], c i double[x], c i int64_t[x], c io dcomplex[xx], c io dcomplex[xx], c io dcomplex[xx], c io dcomplex[xx], c io dcomplex[xx])';
[C1, C2, C3, C4, C5] = AxiStokes3D_mex(mex_id_, ilayer, nt, tx, tnx, nq, sx, snx, M, mu, nrow, C1, C2, C3, C4, C5, 1, 1, nt, nt, 1, nq, nq, 1, 1, 1, nrow, c3m, nrow, c3m, nrow, c3m, nrow, c3m, nrow, c3m);
C1=reshape(C1,nrow,3*nq,Mp1); C2=reshape(C2,nrow,3*nq,Mp1); C3=reshape(C3,nrow,3*nq,Mp1); C4=reshape(C4,nrow,3*nq,Mp1); C5=reshape(C5,nrow,3*nq,Mp1);
end

% ============================================================
% kdtree (axk_): build-once / query-ball, backend = 1 (KD_NBODY) | 2 (KD_TAIYA), required.
%   axk_kdtree_build_mex(backend, pts, leafsize);   % pts is np-by-3 (rows = points, like tin.x')
%   [idx, nout] = axk_kdtree_ball_mex(qpt, radius, nmax);   % idx: 1-based indices within radius
%   axk_kdtree_free_mex();
% ============================================================

% build the tree (state held in the Fortran module until free)
