! Top-level (free-standing) mex-facing wrappers for axissymstok_kernelsplit_mod.
! NOT inside a module.  Prefix axa_ = ax + a (axissymstok module).

! SLP split coefficients.  C1..C5 are (M+1)x3nt x3nq passed as 2D (M+1)x(9 nt nq) from the .mw
! (column-major sequence association: the 3D Fortran view = the 2D mwrap buffer, same memory).
subroutine axa_coef_slp_r64(nt, nq, M, mu, zr, zi, yr, yi, C1, C2, C3, C4, C5)
  use axissymstok_kernelsplit_mod, only: slp => axissymstok_kernelsplit_coef_slp_r64
  implicit none
  integer(8), intent(in)    :: nt, nq, M
  real(8),    intent(in)    :: mu, zr(nt), zi(nt), yr(nq), yi(nq)
  real(8),    intent(inout) :: C1(M+1,3*nt,3*nq), C2(M+1,3*nt,3*nq), C3(M+1,3*nt,3*nq)
  real(8),    intent(inout) :: C4(M+1,3*nt,3*nq), C5(M+1,3*nt,3*nq)
  call slp(nt, nq, M, mu, zr, zi, yr, yi, C1, C2, C3, C4, C5)
end subroutine axa_coef_slp_r64

! DLP split coefficients.  C1,C2 real; C3,C4,C5 complex (carry the normal nu).  Same 2D<->3D
! sequence-association as the SLP wrapper.
subroutine axa_coef_dlp_r64(nt, nq, M, mu, zr, zi, yr, yi, nvr, nvi, C1, C2, C3, C4, C5)
  use axissymstok_kernelsplit_mod, only: dlp => axissymstok_kernelsplit_coef_dlp_r64
  implicit none
  integer(8), intent(in)    :: nt, nq, M
  real(8),    intent(in)    :: mu, zr(nt), zi(nt), yr(nq), yi(nq), nvr(nq), nvi(nq)
  real(8),    intent(inout) :: C1(M+1,3*nt,3*nq), C2(M+1,3*nt,3*nq)
  complex(8), intent(inout) :: C3(M+1,3*nt,3*nq), C4(M+1,3*nt,3*nq), C5(M+1,3*nt,3*nq)
  call dlp(nt, nq, M, mu, zr, zi, yr, yi, nvr, nvi, C1, C2, C3, C4, C5)
end subroutine axa_coef_dlp_r64

! SLPn split coefficients (MER + SWIRL).  TARGET normal nxr,nxi.  MER C1,C2 real, C3,C4 cplx;
! SWIRL (unit n_theta) C1s,C2s real, C3s cplx.
subroutine axa_coef_slpn_r64(nt, nq, M, mu, zr, zi, yr, yi, nxr, nxi, C1, C2, C3, C4, C1s, C2s, C3s)
  use axissymstok_kernelsplit_mod, only: slpn => axissymstok_kernelsplit_coef_slpn_r64
  implicit none
  integer(8), intent(in)    :: nt, nq, M
  real(8),    intent(in)    :: mu, zr(nt), zi(nt), yr(nq), yi(nq), nxr(nt), nxi(nt)
  real(8),    intent(inout) :: C1(M+1,3*nt,3*nq), C2(M+1,3*nt,3*nq)
  complex(8), intent(inout) :: C3(M+1,3*nt,3*nq), C4(M+1,3*nt,3*nq)
  real(8),    intent(inout) :: C1s(M+1,3*nt,3*nq), C2s(M+1,3*nt,3*nq)
  complex(8), intent(inout) :: C3s(M+1,3*nt,3*nq)
  call slpn(nt, nq, M, mu, zr, zi, yr, yi, nxr, nxi, C1, C2, C3, C4, C1s, C2s, C3s)
end subroutine axa_coef_slpn_r64

! Naive modal kernels (no singular quadrature).  K is the plain 3nt x 3ns real block.
subroutine axa_kernel_slp_r64(nt, ns, srcr, srcz, tgtr, tgtz, m, mu, K)
  use axissymstok_kernelsplit_mod, only: kslp => axissymstok_kernel_slp_r64
  implicit none
  integer(8), intent(in)    :: nt, ns, m
  real(8),    intent(in)    :: srcr(ns), srcz(ns), tgtr(nt), tgtz(nt), mu
  real(8),    intent(inout) :: K(3*nt,3*ns)
  call kslp(nt, ns, srcr, srcz, tgtr, tgtz, m, mu, K)
end subroutine axa_kernel_slp_r64

subroutine axa_kernel_dlp_r64(nt, ns, srcr, srcz, srcnr, srcnz, tgtr, tgtz, m, mu, K)
  use axissymstok_kernelsplit_mod, only: kdlp => axissymstok_kernel_dlp_r64
  implicit none
  integer(8), intent(in)    :: nt, ns, m
  real(8),    intent(in)    :: srcr(ns), srcz(ns), srcnr(ns), srcnz(ns), tgtr(nt), tgtz(nt), mu
  real(8),    intent(inout) :: K(3*nt,3*ns)
  call kdlp(nt, ns, srcr, srcz, srcnr, srcnz, tgtr, tgtz, m, mu, K)
end subroutine axa_kernel_dlp_r64

! Naive modal SLPn kernel (mer + swirl); TARGET normal (tgtnr,tgtnz,tgtnth).
subroutine axa_kernel_slpn_r64(nt, ns, srcr, srcz, tgtr, tgtz, tgtnr, tgtnz, tgtnth, m, mu, K, Ksw)
  use axissymstok_kernelsplit_mod, only: kslpn => axissymstok_kernel_slpn_r64
  implicit none
  integer(8), intent(in)    :: nt, ns, m
  real(8),    intent(in)    :: srcr(ns), srcz(ns), tgtr(nt), tgtz(nt), tgtnr(nt), tgtnz(nt), tgtnth(nt), mu
  real(8),    intent(inout) :: K(3*nt,3*ns), Ksw(3*nt,3*ns)
  call kslpn(nt, ns, srcr, srcz, tgtr, tgtz, tgtnr, tgtnz, tgtnth, m, mu, K, Ksw)
end subroutine axa_kernel_slpn_r64
