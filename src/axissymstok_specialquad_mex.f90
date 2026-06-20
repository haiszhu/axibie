! Top-level (free-standing) mex-facing wrappers for axissymstok_specialquad_mod.
! Prefix axss_ = ax + stok specialquad (per NAMING.md).  Complex args pass as dcomplex
! (no real/imag split; mwrap -c99complex), so wrapper + module signatures are identical in
! name and order; each wrapper just forwards to the module compute routine.

subroutine axss_slp_blockmat_r64(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, mu, A)
  use axissymstok_specialquad_mod, only: blk => axissymstok_slp_blockmat_r64
  implicit none
  integer(8), intent(in)    :: nt, p, np, iside, iclosed
  complex(8), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
  real(8),    intent(in)    :: sws(p*np), tpan(np+1), mu
  real(8),    intent(inout) :: A(2*nt, 2*np*p)
  call blk(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, mu, A)
end subroutine axss_slp_blockmat_r64

subroutine axss_slp_blockmat_nmode_r64(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, mu, A)
  use axissymstok_specialquad_mod, only: blk => axissymstok_slp_blockmat_nmode_r64
  implicit none
  integer(8), intent(in)    :: nt, p, np, M, iside, iclosed
  complex(8), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
  real(8),    intent(in)    :: sws(p*np), tpan(np+1), mu
  complex(8), intent(inout) :: A(3*nt, 3*np*p, M+1)
  call blk(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, mu, A)
end subroutine axss_slp_blockmat_nmode_r64

subroutine axss_slpn_blockmat_r64(nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, mu, A)
  use axissymstok_specialquad_mod, only: blk => axissymstok_slpn_blockmat_r64
  implicit none
  integer(8), intent(in)    :: nt, p, np, iside, iclosed
  complex(8), intent(in)    :: tx(nt), tnx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
  real(8),    intent(in)    :: sws(p*np), tpan(np+1), mu
  real(8),    intent(inout) :: A(2*nt, 2*np*p)
  call blk(nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, mu, A)
end subroutine axss_slpn_blockmat_r64

subroutine axss_slpn_blockmat_nmode_r64(nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, mu, A)
  use axissymstok_specialquad_mod, only: blk => axissymstok_slpn_blockmat_nmode_r64
  implicit none
  integer(8), intent(in)    :: nt, p, np, M, iside, iclosed
  complex(8), intent(in)    :: tx(nt), tnx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
  real(8),    intent(in)    :: sws(p*np), tpan(np+1), mu
  complex(8), intent(inout) :: A(3*nt, 3*np*p, M+1)
  call blk(nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, mu, A)
end subroutine axss_slpn_blockmat_nmode_r64

subroutine axss_dlp_blockmat_r64(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, mu, A)
  use axissymstok_specialquad_mod, only: blk => axissymstok_dlp_blockmat_r64
  implicit none
  integer(8), intent(in)    :: nt, p, np, iside, iclosed
  complex(8), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
  real(8),    intent(in)    :: sws(p*np), tpan(np+1), mu
  real(8),    intent(inout) :: A(2*nt, 2*np*p)
  call blk(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, mu, A)
end subroutine axss_dlp_blockmat_r64

subroutine axss_dlp_blockmat_nmode_r64(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, mu, A)
  use axissymstok_specialquad_mod, only: blk => axissymstok_dlp_blockmat_nmode_r64
  implicit none
  integer(8), intent(in)    :: nt, p, np, M, iside, iclosed
  complex(8), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
  real(8),    intent(in)    :: sws(p*np), tpan(np+1), mu
  complex(8), intent(inout) :: A(3*nt, 3*np*p, M+1)
  call blk(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, mu, A)
end subroutine axss_dlp_blockmat_nmode_r64

subroutine axss_dlpn_blockmat_r64(nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, mu, A)
  use axissymstok_specialquad_mod, only: blk => axissymstok_dlpn_blockmat_r64
  implicit none
  integer(8), intent(in)    :: nt, p, np, iside, iclosed
  complex(8), intent(in)    :: tx(nt), tnx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
  real(8),    intent(in)    :: sws(p*np), tpan(np+1), mu
  real(8),    intent(inout) :: A(2*nt, 2*np*p)
  call blk(nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, mu, A)
end subroutine axss_dlpn_blockmat_r64
