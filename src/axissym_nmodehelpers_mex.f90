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
