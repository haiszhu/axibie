! Top-level (free-standing) mex-facing wrappers for axissymlap_kernelsplit_mod.
! NOT inside a module.  Prefix axl_ = ax + l (axissymlap kernelsplit; 0th azimuthal mode).
!
! All split coefficients are REAL nt x nq densities (the normals nu=t.nx, nu'=s.nx live in
! the special-quad matrices, not the coefficients).  Naming follows the Stokes axa_coef_*
! wrappers with the 0th-mode "0" marker.  Thin delegation to the secret module procs.

! ---- SLP  S : C1 log, C2 smooth ----
subroutine axl_coef_slp0_r64(nt, nq, zr, zi, yr, yi, C1, C2)
  use axissymlap_kernelsplit_mod, only: slp0 => axissymlap_kernelsplit_coef_slp0_r64
  implicit none
  integer(8), intent(in)    :: nt, nq
  real(8),    intent(in)    :: zr(nt), zi(nt), yr(nq), yi(nq)
  real(8),    intent(inout) :: C1(nt,nq), C2(nt,nq)
  call slp0(nt, nq, zr, zi, yr, yi, C1, C2)
end subroutine axl_coef_slp0_r64

! ---- SLPn  S' : target normal (nxr,nxi).  C1 log, C2 smooth, C3 target-Cauchy density ----
subroutine axl_coef_slpn0_r64(nt, nq, zr, zi, yr, yi, nxr, nxi, C1, C2, C3)
  use axissymlap_kernelsplit_mod, only: slpn0 => axissymlap_kernelsplit_coef_slpn0_r64
  implicit none
  integer(8), intent(in)    :: nt, nq
  real(8),    intent(in)    :: zr(nt), zi(nt), yr(nq), yi(nq), nxr(nt), nxi(nt)
  real(8),    intent(inout) :: C1(nt,nq), C2(nt,nq), C3(nt,nq)
  call slpn0(nt, nq, zr, zi, yr, yi, nxr, nxi, C1, C2, C3)
end subroutine axl_coef_slpn0_r64

! ---- DLP  D : source normal (nvr,nvi).  C1 log, C2 smooth, C3 source-Cauchy density ----
subroutine axl_coef_dlp0_r64(nt, nq, zr, zi, yr, yi, nvr, nvi, C1, C2, C3)
  use axissymlap_kernelsplit_mod, only: dlp0 => axissymlap_kernelsplit_coef_dlp0_r64
  implicit none
  integer(8), intent(in)    :: nt, nq
  real(8),    intent(in)    :: zr(nt), zi(nt), yr(nq), yi(nq), nvr(nq), nvi(nq)
  real(8),    intent(inout) :: C1(nt,nq), C2(nt,nq), C3(nt,nq)
  call dlp0(nt, nq, zr, zi, yr, yi, nvr, nvi, C1, C2, C3)
end subroutine axl_coef_dlp0_r64

! ---- DLPn  D' : both normals.  C1 log, C2 smooth, C3a target-Cauchy, C3b source-Cauchy,
!      C4 (nu nu') deriv-Cauchy density ----
subroutine axl_coef_dlpn0_r64(nt, nq, zr, zi, yr, yi, nxr, nxi, nvr, nvi, C1, C2, C3a, C3b, C4)
  use axissymlap_kernelsplit_mod, only: dlpn0 => axissymlap_kernelsplit_coef_dlpn0_r64
  implicit none
  integer(8), intent(in)    :: nt, nq
  real(8),    intent(in)    :: zr(nt), zi(nt), yr(nq), yi(nq), nxr(nt), nxi(nt), nvr(nq), nvi(nq)
  real(8),    intent(inout) :: C1(nt,nq), C2(nt,nq), C3a(nt,nq), C3b(nt,nq), C4(nt,nq)
  call dlpn0(nt, nq, zr, zi, yr, yi, nxr, nxi, nvr, nvi, C1, C2, C3a, C3b, C4)
end subroutine axl_coef_dlpn0_r64

! ---- DLPnn  D'' : both normals.  C1 log, C2 smooth, C3a target-Cauchy, C3b source-Cauchy,
!      C4a (nu^2) deriv-Cauchy, C4b (nu nu') deriv-Cauchy, C5 cubic density ----
subroutine axl_coef_dlpnn0_r64(nt, nq, zr, zi, yr, yi, nxr, nxi, nvr, nvi, C1, C2, C3a, C3b, C4a, C4b, C5)
  use axissymlap_kernelsplit_mod, only: dlpnn0 => axissymlap_kernelsplit_coef_dlpnn0_r64
  implicit none
  integer(8), intent(in)    :: nt, nq
  real(8),    intent(in)    :: zr(nt), zi(nt), yr(nq), yi(nq), nxr(nt), nxi(nt), nvr(nq), nvi(nq)
  real(8),    intent(inout) :: C1(nt,nq), C2(nt,nq), C3a(nt,nq), C3b(nt,nq), C4a(nt,nq), C4b(nt,nq), C5(nt,nq)
  call dlpnn0(nt, nq, zr, zi, yr, yi, nxr, nxi, nvr, nvi, C1, C2, C3a, C3b, C4a, C4b, C5)
end subroutine axl_coef_dlpnn0_r64
