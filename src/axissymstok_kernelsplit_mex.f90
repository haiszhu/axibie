! Top-level mex-facing wrappers for axissymstok_kernelsplit_mod (prefix axsk_).
! The n-mode close-eval split coefficients, exposed so the pure-MATLAB n-mode block builders
! (utils/axss_*_blockmat_nmode.m) can call the same compiled split as the Fortran builders.
! Each wrapper just forwards to the (secret) module compute routine.

subroutine axsk_slp_coef_nmode_r64(nt, tx, nq, sx, M, mu, C1, C2, C3)
  use axissymstok_kernelsplit_mod, only: cf => axissymstok_slp_coef_nmode_r64
  implicit none
  integer(8), intent(in)    :: nt, nq, M
  complex(8), intent(in)    :: tx(nt), sx(nq)
  real(8),    intent(in)    :: mu
  complex(8), intent(inout) :: C1(3*nt,3*nq,M+1), C2(3*nt,3*nq,M+1), C3(3*nt,3*nq,M+1)
  call cf(nt, tx, nq, sx, M, mu, C1, C2, C3)
end subroutine axsk_slp_coef_nmode_r64

subroutine axsk_slpn_coef_nmode_r64(nt, tx, tnx, nq, sx, M, mu, C1, C2, C3, C4)
  use axissymstok_kernelsplit_mod, only: cf => axissymstok_slpn_coef_nmode_r64
  implicit none
  integer(8), intent(in)    :: nt, nq, M
  complex(8), intent(in)    :: tx(nt), tnx(nt), sx(nq)
  real(8),    intent(in)    :: mu
  complex(8), intent(inout) :: C1(3*nt,3*nq,M+1), C2(3*nt,3*nq,M+1), C3(3*nt,3*nq,M+1), C4(3*nt,3*nq,M+1)
  call cf(nt, tx, tnx, nq, sx, M, mu, C1, C2, C3, C4)
end subroutine axsk_slpn_coef_nmode_r64

subroutine axsk_dlp_coef_nmode_r64(nt, tx, nq, sx, snx, M, mu, C1, C2, C3, C4)
  use axissymstok_kernelsplit_mod, only: cf => axissymstok_dlp_coef_nmode_r64
  implicit none
  integer(8), intent(in)    :: nt, nq, M
  complex(8), intent(in)    :: tx(nt), sx(nq), snx(nq)
  real(8),    intent(in)    :: mu
  complex(8), intent(inout) :: C1(3*nt,3*nq,M+1), C2(3*nt,3*nq,M+1), C3(3*nt,3*nq,M+1), C4(3*nt,3*nq,M+1)
  call cf(nt, tx, nq, sx, snx, M, mu, C1, C2, C3, C4)
end subroutine axsk_dlp_coef_nmode_r64

subroutine axsk_dlpn_coef_nmode_r64(nt, tx, tnx, nq, sx, snx, M, mu, C1, C2, C3, C4, C5)
  use axissymstok_kernelsplit_mod, only: cf => axissymstok_dlpn_coef_nmode_r64
  implicit none
  integer(8), intent(in)    :: nt, nq, M
  complex(8), intent(in)    :: tx(nt), tnx(nt), sx(nq), snx(nq)
  real(8),    intent(in)    :: mu
  complex(8), intent(inout) :: C1(3*nt,3*nq,M+1), C2(3*nt,3*nq,M+1), C3(3*nt,3*nq,M+1), C4(3*nt,3*nq,M+1), C5(3*nt,3*nq,M+1)
  call cf(nt, tx, tnx, nq, sx, snx, M, mu, C1, C2, C3, C4, C5)
end subroutine axsk_dlpn_coef_nmode_r64

subroutine axsk_slppres_coef_nmode_r64(nt, tx, nq, sx, M, C1, C2, C3)
  use axissymstok_kernelsplit_mod, only: cf => axissymstok_slppres_coef_nmode_r64
  implicit none
  integer(8), intent(in)    :: nt, nq, M
  complex(8), intent(in)    :: tx(nt), sx(nq)
  complex(8), intent(inout) :: C1(nt,3*nq,M+1), C2(nt,3*nq,M+1), C3(nt,3*nq,M+1)
  call cf(nt, tx, nq, sx, M, C1, C2, C3)
end subroutine axsk_slppres_coef_nmode_r64

subroutine axsk_dlppres_coef_nmode_r64(nt, tx, nq, sx, snx, M, C1, C2, C3, C4)
  use axissymstok_kernelsplit_mod, only: cf => axissymstok_dlppres_coef_nmode_r64
  implicit none
  integer(8), intent(in)    :: nt, nq, M
  complex(8), intent(in)    :: tx(nt), sx(nq), snx(nq)
  complex(8), intent(inout) :: C1(nt,3*nq,M+1), C2(nt,3*nq,M+1), C3(nt,3*nq,M+1), C4(nt,3*nq,M+1)
  call cf(nt, tx, nq, sx, snx, M, C1, C2, C3, C4)
end subroutine axsk_dlppres_coef_nmode_r64

subroutine axsk_coef_r64(ilayer, nt, tx, tnx, nq, sx, snx, M, mu, nrow, C1, C2, C3, C4, C5)
  use axissymstok_kernelsplit_mod, only: cf => axissymstok_coef_r64
  implicit none
  integer(8), intent(in)    :: ilayer, nt, nq, M, nrow
  complex(8), intent(in)    :: tx(nt), tnx(nt), sx(nq), snx(nq)
  real(8),    intent(in)    :: mu
  complex(8), intent(inout) :: C1(nrow,3*nq,M+1), C2(nrow,3*nq,M+1), C3(nrow,3*nq,M+1), C4(nrow,3*nq,M+1), C5(nrow,3*nq,M+1)
  call cf(ilayer, nt, tx, tnx, nq, sx, snx, M, mu, nrow, C1, C2, C3, C4, C5)
end subroutine axsk_coef_r64
