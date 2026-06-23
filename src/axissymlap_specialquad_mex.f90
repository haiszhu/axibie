! Top-level (free-standing) mex-facing wrappers for axissymlap_specialquad_mod.
! Prefix axls_ = ax + lap specialquad (per NAMING.md).  Complex args pass as dcomplex
! (no real/imag split; mwrap -c99complex), so wrapper + module signatures are identical in
! name and order; each wrapper just forwards to the module compute routine.  Laplace is
! scalar (no mu); the n-mode block A is REAL (no phi-cross imaginary part).

subroutine axls_slp_blockmat_r64(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, A)
  use axissymlap_specialquad_mod, only: blk => axissymlap_slp_blockmat_r64
  implicit none
  integer(8), intent(in)    :: nt, p, np, iside, iclosed
  complex(8), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
  real(8),    intent(in)    :: sws(p*np), tpan(np+1)
  real(8),    intent(inout) :: A(nt, np*p)
  call blk(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, A)
end subroutine axls_slp_blockmat_r64

subroutine axls_slp_blockmat_nmode_r64(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, A)
  use axissymlap_specialquad_mod, only: blk => axissymlap_slp_blockmat_nmode_r64
  implicit none
  integer(8), intent(in)    :: nt, p, np, M, iside, iclosed
  complex(8), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
  real(8),    intent(in)    :: sws(p*np), tpan(np+1)
  real(8),    intent(inout) :: A(nt, np*p, M+1)
  call blk(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, A)
end subroutine axls_slp_blockmat_nmode_r64

subroutine axls_slpn_blockmat_r64(nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, A)
  use axissymlap_specialquad_mod, only: blk => axissymlap_slpn_blockmat_r64
  implicit none
  integer(8), intent(in)    :: nt, p, np, iside, iclosed
  complex(8), intent(in)    :: tx(nt), tnx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
  real(8),    intent(in)    :: sws(p*np), tpan(np+1)
  real(8),    intent(inout) :: A(nt, np*p)
  call blk(nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, A)
end subroutine axls_slpn_blockmat_r64

subroutine axls_slpn_blockmat_nmode_r64(nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, A)
  use axissymlap_specialquad_mod, only: blk => axissymlap_slpn_blockmat_nmode_r64
  implicit none
  integer(8), intent(in)    :: nt, p, np, M, iside, iclosed
  complex(8), intent(in)    :: tx(nt), tnx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
  real(8),    intent(in)    :: sws(p*np), tpan(np+1)
  real(8),    intent(inout) :: A(nt, np*p, M+1)
  call blk(nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, A)
end subroutine axls_slpn_blockmat_nmode_r64

subroutine axls_dlp_blockmat_r64(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, A)
  use axissymlap_specialquad_mod, only: blk => axissymlap_dlp_blockmat_r64
  implicit none
  integer(8), intent(in)    :: nt, p, np, iside, iclosed
  complex(8), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
  real(8),    intent(in)    :: sws(p*np), tpan(np+1)
  real(8),    intent(inout) :: A(nt, np*p)
  call blk(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, A)
end subroutine axls_dlp_blockmat_r64

subroutine axls_dlp_blockmat_nmode_r64(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, A)
  use axissymlap_specialquad_mod, only: blk => axissymlap_dlp_blockmat_nmode_r64
  implicit none
  integer(8), intent(in)    :: nt, p, np, M, iside, iclosed
  complex(8), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
  real(8),    intent(in)    :: sws(p*np), tpan(np+1)
  real(8),    intent(inout) :: A(nt, np*p, M+1)
  call blk(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, A)
end subroutine axls_dlp_blockmat_nmode_r64

subroutine axls_dlpn_blockmat_r64(nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, A)
  use axissymlap_specialquad_mod, only: blk => axissymlap_dlpn_blockmat_r64
  implicit none
  integer(8), intent(in)    :: nt, p, np, iside, iclosed
  complex(8), intent(in)    :: tx(nt), tnx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
  real(8),    intent(in)    :: sws(p*np), tpan(np+1)
  real(8),    intent(inout) :: A(nt, np*p)
  call blk(nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, A)
end subroutine axls_dlpn_blockmat_r64

subroutine axls_dlpn_blockmat_nmode_r64(nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, A)
  use axissymlap_specialquad_mod, only: blk => axissymlap_dlpn_blockmat_nmode_r64
  implicit none
  integer(8), intent(in)    :: nt, p, np, M, iside, iclosed
  complex(8), intent(in)    :: tx(nt), tnx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
  real(8),    intent(in)    :: sws(p*np), tpan(np+1)
  real(8),    intent(inout) :: A(nt, np*p, M+1)
  call blk(nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, A)
end subroutine axls_dlpn_blockmat_nmode_r64
