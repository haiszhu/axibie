function [C1, C2, C3, C4] = axsk_dlppres_coef_nmode_mex(nt, tx, nq, sx, snx, M)
nt=double(nt); nq=double(nq); M=double(M); tx=double(tx(:)); sx=double(sx(:)); snx=double(snx(:));
Mp1=M+1; c3m=3*nq*Mp1;
C1=complex(zeros(nt,c3m)); C2=complex(zeros(nt,c3m)); C3=complex(zeros(nt,c3m)); C4=complex(zeros(nt,c3m));
mex_id_ = 'axsk_dlppres_coef_nmode_r64(c i int64_t[x], c i dcomplex[x], c i int64_t[x], c i dcomplex[x], c i dcomplex[x], c i int64_t[x], c io dcomplex[xx], c io dcomplex[xx], c io dcomplex[xx], c io dcomplex[xx])';
[C1, C2, C3, C4] = AxiStokes3D_mex(mex_id_, nt, tx, nq, sx, snx, M, C1, C2, C3, C4, 1, nt, 1, nq, nq, 1, nt, c3m, nt, c3m, nt, c3m, nt, c3m);
C1=reshape(C1,nt,3*nq,Mp1); C2=reshape(C2,nt,3*nq,Mp1); C3=reshape(C3,nt,3*nq,Mp1); C4=reshape(C4,nt,3*nq,Mp1);
end

% LEVEL-1 MASTER (Stokes): coef dispatch, ilayer 1 SLP / 2 SLPn / 3 DLP / 4 DLPn / 5 SLPpres / 6 DLPpres.
% Buckets (nrow, 3*nq, M+1): nrow = 3*nt for ilayer 1-4, nt for the pressure rows 5-6.  Unused slots zeroed.
