! Top-level mex-facing wrappers for the numerics helpers the n-mode block builders call:
! Gauss-Legendre + spectral diff (gauss), Lagrange interpolation matrix (lagrange_interp),
! the all-modes Q-side carrier (modal_green_all_far), and Helsing close-eval special quad
! (sdspecialquad).  Exposed so the pure-MATLAB n-mode builders (utils/axls_*_blockmat_nmode.m)
! call the SAME compiled routines as the Fortran builders.  Each wrapper just forwards.

subroutine axt_gauss_r64(n, tgl, wgl, Dgl)
  use axistokes3d_mod, only: g => gauss_r64
  implicit none
  integer(8), intent(in)    :: n
  real(8),    intent(inout) :: tgl(n), wgl(n), Dgl(n,n)
  call g(n, tgl, wgl, Dgl)
end subroutine axt_gauss_r64

subroutine axt_lagrange_interp_r64(ns, xs, nt, xt, Mat)
  use axistokes3d_mod, only: li => lagrange_interp_r64
  implicit none
  integer(8), intent(in)    :: ns, nt
  real(8),    intent(in)    :: xs(ns), xt(nt)
  real(8),    intent(inout) :: Mat(nt,ns)
  call li(ns, xs, nt, xt, Mat)
end subroutine axt_lagrange_interp_r64

subroutine axt_lagrange_eval_r64(ns, xs, yc, nt, xt, yout)
  ! Evaluate Lagrange interpolant: yout(j) = sum_k Mat(j,k)*yc(k), j=1..nt
  ! where Mat = lagrange_interp(ns,xs,nt,xt).  matmul done in Fortran so the
  ! arithmetic matches the block builders (same as their internal ec2=matmul(IPe2,xc)).
  use axistokes3d_mod, only: li => lagrange_interp_r64
  implicit none
  integer(8),  intent(in)    :: ns, nt
  real(8),     intent(in)    :: xs(ns), xt(nt)
  complex(8),  intent(in)    :: yc(ns)
  complex(8),  intent(inout) :: yout(nt)
  real(8) :: Mat(nt,ns)
  call li(ns, xs, nt, xt, Mat)
  yout = matmul(Mat, yc)
end subroutine axt_lagrange_eval_r64

subroutine axt_rmatcz_r64(m, n, A, x, y)
  ! y = matmul(A, x), A real(m,n), x complex(n) -> y complex(m).
  ! Fortran matmul with -ffp-contract=off, matching the block builders' arithmetic.
  implicit none
  integer(8),  intent(in)    :: m, n
  real(8),     intent(in)    :: A(m,n)
  complex(8),  intent(in)    :: x(n)
  complex(8),  intent(inout) :: y(m)
  y = matmul(A, x)
end subroutine axt_rmatcz_r64

subroutine axmg_modal_green_all_far_r64(chi, M, vk, ve)
  use axisym_modal_green_mod, only: mg => modal_green_all_far_r64
  implicit none
  integer(8), intent(in)    :: M
  real(8),    intent(in)    :: chi
  real(8),    intent(inout) :: vk(M+1), ve(M+1)
  call mg(chi, M, vk, ve)
end subroutine axmg_modal_green_all_far_r64

subroutine axsq_sdspecialquad_r64(nt, zt, p, zsrc, nzsrc, wzp, za, zb, iside, As, Ad, A1, A2, A3, A4)
  use specialquad_mod, only: sq => sdspecialquad_r64
  implicit none
  integer(8), intent(in)    :: nt, p, iside
  complex(8), intent(in)    :: zt(nt), zsrc(p), nzsrc(p), wzp(p), za, zb
  real(8),    intent(inout) :: As(p,nt), A1(p,nt), A2(p,nt), A3(p,nt), A4(p,nt)
  complex(8), intent(inout) :: Ad(p,nt)
  call sq(nt, zt, p, zsrc, nzsrc, wzp, za, zb, iside, As, Ad, A1, A2, A3, A4)
end subroutine axsq_sdspecialquad_r64

subroutine axt_sto3dslp_eval_r64(nt, tx, ns, sx, sw, f, u)
  use axistokes3d_mod, only: ev => sto3dslp_eval_r64
  implicit none
  integer(8), intent(in)    :: nt, ns
  real(8),    intent(in)    :: tx(3,nt), sx(3,ns), sw(ns), f(ns,3)
  real(8),    intent(inout) :: u(3,nt)
  call ev(nt, tx, ns, sx, sw, f, u)
end subroutine axt_sto3dslp_eval_r64

subroutine axt_sto3ddlp_eval_r64(nt, tx, ns, sx, snx, sw, f, u)
  use axistokes3d_mod, only: ev => sto3ddlp_eval_r64
  implicit none
  integer(8), intent(in)    :: nt, ns
  real(8),    intent(in)    :: tx(3,nt), sx(3,ns), snx(3,ns), sw(ns), f(ns,3)
  real(8),    intent(inout) :: u(3,nt)
  call ev(nt, tx, ns, sx, snx, sw, f, u)
end subroutine axt_sto3ddlp_eval_r64
