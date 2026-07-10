! Top-level mex-facing wrappers for axissymlap_kernelsplit_mod (prefix axlk_).
! The n-mode close-eval split coefficients, exposed so the pure-MATLAB n-mode block builders
! (utils/axls_*_blockmat_nmode.m) can call the same compiled split as the Fortran builders.
! Each wrapper just forwards to the (secret) module compute routine.

subroutine axlk_slp_coef_nmode_r64(nt, tx, nq, sx, M, C1, C2)
  use axissymlap_kernelsplit_mod, only: cf => axissymlap_slp_coef_nmode_r64
  implicit none
  integer(8), intent(in)    :: nt, nq, M
  complex(8), intent(in)    :: tx(nt), sx(nq)
  real(8),    intent(inout) :: C1(nt,nq,M+1), C2(nt,nq,M+1)
  call cf(nt, tx, nq, sx, M, C1, C2)
end subroutine axlk_slp_coef_nmode_r64

subroutine axlk_slpn_coef_nmode_r64(nt, tx, tnx, nq, sx, M, C1, C2, C3)
  use axissymlap_kernelsplit_mod, only: cf => axissymlap_slpn_coef_nmode_r64
  implicit none
  integer(8), intent(in)    :: nt, nq, M
  complex(8), intent(in)    :: tx(nt), tnx(nt), sx(nq)
  real(8),    intent(inout) :: C1(nt,nq,M+1), C2(nt,nq,M+1), C3(nt,nq,M+1)
  call cf(nt, tx, tnx, nq, sx, M, C1, C2, C3)
end subroutine axlk_slpn_coef_nmode_r64

subroutine axlk_dlp_coef_nmode_r64(nt, tx, nq, sx, snx, M, C1, C2, C3)
  use axissymlap_kernelsplit_mod, only: cf => axissymlap_dlp_coef_nmode_r64
  implicit none
  integer(8), intent(in)    :: nt, nq, M
  complex(8), intent(in)    :: tx(nt), sx(nq), snx(nq)
  real(8),    intent(inout) :: C1(nt,nq,M+1), C2(nt,nq,M+1), C3(nt,nq,M+1)
  call cf(nt, tx, nq, sx, snx, M, C1, C2, C3)
end subroutine axlk_dlp_coef_nmode_r64

subroutine axlk_dlpn_coef_nmode_r64(nt, tx, tnx, nq, sx, snx, M, C1, C2, C3a, C3b, C4)
  use axissymlap_kernelsplit_mod, only: cf => axissymlap_dlpn_coef_nmode_r64
  implicit none
  integer(8), intent(in)    :: nt, nq, M
  complex(8), intent(in)    :: tx(nt), tnx(nt), sx(nq), snx(nq)
  real(8),    intent(inout) :: C1(nt,nq,M+1), C2(nt,nq,M+1), C3a(nt,nq,M+1), C3b(nt,nq,M+1), C4(nt,nq,M+1)
  call cf(nt, tx, tnx, nq, sx, snx, M, C1, C2, C3a, C3b, C4)
end subroutine axlk_dlpn_coef_nmode_r64

subroutine axlk_coef_r64(ilayer, nt, tx, tnx, nq, sx, snx, M, C1, C2, C3a, C3b, C4)
  use axissymlap_kernelsplit_mod, only: cf => axissymlap_coef_r64
  implicit none
  integer(8), intent(in)    :: ilayer, nt, nq, M
  complex(8), intent(in)    :: tx(nt), tnx(nt), sx(nq), snx(nq)
  real(8),    intent(inout) :: C1(nt,nq,M+1), C2(nt,nq,M+1), C3a(nt,nq,M+1), C3b(nt,nq,M+1), C4(nt,nq,M+1)
  call cf(ilayer, nt, tx, tnx, nq, sx, snx, M, C1, C2, C3a, C3b, C4)
end subroutine axlk_coef_r64
