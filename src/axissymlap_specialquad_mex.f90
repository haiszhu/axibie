! Top-level (free-standing) mex-facing wrappers for axissymlap_specialquad_mod.
! Prefix axls_ = ax + ls (axissymlap specialquad; 0th azimuthal mode close-eval operator).
! Complex inputs split into real/imag (same convention as axm_); W is the nt x p REAL block.
! Scalar (no 3x interleave), 0th-mode (no M), Laplace (no mu).

subroutine axls_specialquad_slp0_r64(nt, p, hm, iside, zr, zi, yr, yi, ypr, ypi, &
    zar, zai, zbr, zbi, gw, Pmat, W)
  use axissymlap_specialquad_mod, only: slp => axissymlap_specialquad_slp0_r64
  implicit none
  integer(8), intent(in)    :: nt, p, iside
  real(8),    intent(in)    :: hm, zar, zai, zbr, zbi
  real(8),    intent(in)    :: zr(nt), zi(nt), yr(2*p), yi(2*p), ypr(2*p), ypi(2*p)
  real(8),    intent(in)    :: gw(2*p), Pmat(2*p,p)
  real(8),    intent(inout) :: W(nt,p)
  complex(8) :: z(nt), y(2*p), yp(2*p), za, zb
  integer(8) :: k
  do k = 1, nt;   z(k)  = cmplx(zr(k),  zi(k),  8); end do
  do k = 1, 2*p;  y(k)  = cmplx(yr(k),  yi(k),  8); yp(k) = cmplx(ypr(k), ypi(k), 8); end do
  za = cmplx(zar, zai, 8); zb = cmplx(zbr, zbi, 8)
  call slp(nt, z, p, hm, y, yp, za, zb, gw, Pmat, iside, W)
end subroutine axls_specialquad_slp0_r64

subroutine axls_specialquad_slpn0_r64(nt, p, hm, iside, zr, zi, nxr, nxi, yr, yi, ypr, ypi, &
    zar, zai, zbr, zbi, gw, Pmat, W)
  use axissymlap_specialquad_mod, only: slpn => axissymlap_specialquad_slpn0_r64
  implicit none
  integer(8), intent(in)    :: nt, p, iside
  real(8),    intent(in)    :: hm, zar, zai, zbr, zbi
  real(8),    intent(in)    :: zr(nt), zi(nt), nxr(nt), nxi(nt)
  real(8),    intent(in)    :: yr(2*p), yi(2*p), ypr(2*p), ypi(2*p), gw(2*p), Pmat(2*p,p)
  real(8),    intent(inout) :: W(nt,p)
  complex(8) :: z(nt), nx(nt), y(2*p), yp(2*p), za, zb
  integer(8) :: k
  do k = 1, nt;   z(k)  = cmplx(zr(k),  zi(k),  8); nx(k) = cmplx(nxr(k), nxi(k), 8); end do
  do k = 1, 2*p;  y(k)  = cmplx(yr(k),  yi(k),  8); yp(k) = cmplx(ypr(k), ypi(k), 8); end do
  za = cmplx(zar, zai, 8); zb = cmplx(zbr, zbi, 8)
  call slpn(nt, z, nx, p, hm, y, yp, za, zb, gw, Pmat, iside, W)
end subroutine axls_specialquad_slpn0_r64

subroutine axls_specialquad_dlp0_r64(nt, p, hm, iside, zr, zi, yr, yi, ypr, ypi, &
    zar, zai, zbr, zbi, gw, Pmat, W)
  use axissymlap_specialquad_mod, only: dlp => axissymlap_specialquad_dlp0_r64
  implicit none
  integer(8), intent(in)    :: nt, p, iside
  real(8),    intent(in)    :: hm, zar, zai, zbr, zbi
  real(8),    intent(in)    :: zr(nt), zi(nt), yr(2*p), yi(2*p), ypr(2*p), ypi(2*p)
  real(8),    intent(in)    :: gw(2*p), Pmat(2*p,p)
  real(8),    intent(inout) :: W(nt,p)
  complex(8) :: z(nt), y(2*p), yp(2*p), za, zb
  integer(8) :: k
  do k = 1, nt;   z(k)  = cmplx(zr(k),  zi(k),  8); end do
  do k = 1, 2*p;  y(k)  = cmplx(yr(k),  yi(k),  8); yp(k) = cmplx(ypr(k), ypi(k), 8); end do
  za = cmplx(zar, zai, 8); zb = cmplx(zbr, zbi, 8)
  call dlp(nt, z, p, hm, y, yp, za, zb, gw, Pmat, iside, W)
end subroutine axls_specialquad_dlp0_r64

subroutine axls_specialquad_dlpn0_r64(nt, p, hm, iside, zr, zi, nxr, nxi, yr, yi, ypr, ypi, &
    zar, zai, zbr, zbi, gw, Pmat, W)
  use axissymlap_specialquad_mod, only: dlpn => axissymlap_specialquad_dlpn0_r64
  implicit none
  integer(8), intent(in)    :: nt, p, iside
  real(8),    intent(in)    :: hm, zar, zai, zbr, zbi
  real(8),    intent(in)    :: zr(nt), zi(nt), nxr(nt), nxi(nt)
  real(8),    intent(in)    :: yr(2*p), yi(2*p), ypr(2*p), ypi(2*p), gw(2*p), Pmat(2*p,p)
  real(8),    intent(inout) :: W(nt,p)
  complex(8) :: z(nt), nx(nt), y(2*p), yp(2*p), za, zb
  integer(8) :: k
  do k = 1, nt;   z(k)  = cmplx(zr(k),  zi(k),  8); nx(k) = cmplx(nxr(k), nxi(k), 8); end do
  do k = 1, 2*p;  y(k)  = cmplx(yr(k),  yi(k),  8); yp(k) = cmplx(ypr(k), ypi(k), 8); end do
  za = cmplx(zar, zai, 8); zb = cmplx(zbr, zbi, 8)
  call dlpn(nt, z, nx, p, hm, y, yp, za, zb, gw, Pmat, iside, W)
end subroutine axls_specialquad_dlpn0_r64

subroutine axls_specialquad_dlpnn0_r64(nt, p, hm, iside, zr, zi, nxr, nxi, yr, yi, ypr, ypi, &
    zar, zai, zbr, zbi, gw, Pmat, W)
  use axissymlap_specialquad_mod, only: dlpnn => axissymlap_specialquad_dlpnn0_r64
  implicit none
  integer(8), intent(in)    :: nt, p, iside
  real(8),    intent(in)    :: hm, zar, zai, zbr, zbi
  real(8),    intent(in)    :: zr(nt), zi(nt), nxr(nt), nxi(nt)
  real(8),    intent(in)    :: yr(2*p), yi(2*p), ypr(2*p), ypi(2*p), gw(2*p), Pmat(2*p,p)
  real(8),    intent(inout) :: W(nt,p)
  complex(8) :: z(nt), nx(nt), y(2*p), yp(2*p), za, zb
  integer(8) :: k
  do k = 1, nt;   z(k)  = cmplx(zr(k),  zi(k),  8); nx(k) = cmplx(nxr(k), nxi(k), 8); end do
  do k = 1, 2*p;  y(k)  = cmplx(yr(k),  yi(k),  8); yp(k) = cmplx(ypr(k), ypi(k), 8); end do
  za = cmplx(zar, zai, 8); zb = cmplx(zbr, zbi, 8)
  call dlpnn(nt, z, nx, p, hm, y, yp, za, zb, gw, Pmat, iside, W)
end subroutine axls_specialquad_dlpnn0_r64
