! Top-level (free-standing) mex-facing wrappers for axissymstok_specialquad_mod.
! Prefix axm_ = ax + m (modal close-eval operator).
! Complex inputs are passed split into real/imag (same convention as the coef wrappers);
! W is (M+1)x3nt x3p passed as 2D (M+1)x(9 nt p) via column-major sequence association.

subroutine axm_specialquad_dlp_r64(nt, p, M, mu, hm, iside, zr, zi, yr, yi, ypr, ypi, &
    zar, zai, zbr, zbi, gw, Pmat, W)
  use axissymstok_specialquad_mod, only: dlp => axissymstok_specialquad_dlp_r64
  implicit none
  integer(8), intent(in)    :: nt, p, M, iside
  real(8),    intent(in)    :: mu, hm, zar, zai, zbr, zbi
  real(8),    intent(in)    :: zr(nt), zi(nt), yr(2*p), yi(2*p), ypr(2*p), ypi(2*p)
  real(8),    intent(in)    :: gw(2*p), Pmat(2*p,p)
  real(8),    intent(inout) :: W(M+1,3*nt,3*p)
  complex(8) :: z(nt), y(2*p), yp(2*p), za, zb
  integer(8) :: k
  do k = 1, nt
    z(k) = cmplx(zr(k), zi(k), 8)
  end do
  do k = 1, 2*p
    y(k)  = cmplx(yr(k),  yi(k),  8)
    yp(k) = cmplx(ypr(k), ypi(k), 8)
  end do
  za = cmplx(zar, zai, 8); zb = cmplx(zbr, zbi, 8)
  call dlp(nt, z, p, hm, y, yp, za, zb, gw, Pmat, M, mu, iside, W)
end subroutine axm_specialquad_dlp_r64

subroutine axm_specialquad_slp_r64(nt, p, M, mu, hm, iside, zr, zi, yr, yi, ypr, ypi, &
    zar, zai, zbr, zbi, gw, Pmat, W)
  use axissymstok_specialquad_mod, only: slp => axissymstok_specialquad_slp_r64
  implicit none
  integer(8), intent(in)    :: nt, p, M, iside
  real(8),    intent(in)    :: mu, hm, zar, zai, zbr, zbi
  real(8),    intent(in)    :: zr(nt), zi(nt), yr(2*p), yi(2*p), ypr(2*p), ypi(2*p)
  real(8),    intent(in)    :: gw(2*p), Pmat(2*p,p)
  real(8),    intent(inout) :: W(M+1,3*nt,3*p)
  complex(8) :: z(nt), y(2*p), yp(2*p), za, zb
  integer(8) :: k
  do k = 1, nt
    z(k) = cmplx(zr(k), zi(k), 8)
  end do
  do k = 1, 2*p
    y(k)  = cmplx(yr(k),  yi(k),  8)
    yp(k) = cmplx(ypr(k), ypi(k), 8)
  end do
  za = cmplx(zar, zai, 8); zb = cmplx(zbr, zbi, 8)
  call slp(nt, z, p, hm, y, yp, za, zb, gw, Pmat, M, mu, iside, W)
end subroutine axm_specialquad_slp_r64
