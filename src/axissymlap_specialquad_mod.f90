module axissymlap_specialquad_mod
  use axilaplace3d_mod, only: r64
  implicit none
  private
  public :: axissymlap_slp_blockmat_r64
  public :: axissymlap_slp_blockmat_nmode_r64
  public :: axissymlap_slpn_blockmat_r64
  public :: axissymlap_slpn_blockmat_nmode_r64
  public :: axissymlap_dlp_blockmat_r64
  public :: axissymlap_dlp_blockmat_nmode_r64
  public :: axissymlap_dlpn_blockmat_r64
  public :: axissymlap_dlpn_blockmat_nmode_r64
contains

  subroutine axissymlap_slp_blockmat_r64(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, A)
    integer(8),   intent(in)    :: nt, p, np, iside, iclosed
    complex(r64), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: sws(p*np), tpan(np+1)
    real(r64),    intent(inout) :: A(nt, np*p)
  end subroutine axissymlap_slp_blockmat_r64

  subroutine axissymlap_slp_blockmat_nmode_r64(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, A)
    integer(8),   intent(in)    :: nt, p, np, M, iside, iclosed
    complex(r64), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: sws(p*np), tpan(np+1)
    real(r64),    intent(inout) :: A(nt, np*p, M+1)
  end subroutine axissymlap_slp_blockmat_nmode_r64

  subroutine axissymlap_slpn_blockmat_r64(nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, A)
    integer(8),   intent(in)    :: nt, p, np, iside, iclosed
    complex(r64), intent(in)    :: tx(nt), tnx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: sws(p*np), tpan(np+1)
    real(r64),    intent(inout) :: A(nt, np*p)
  end subroutine axissymlap_slpn_blockmat_r64

  subroutine axissymlap_slpn_blockmat_nmode_r64(nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, A)
    integer(8),   intent(in)    :: nt, p, np, M, iside, iclosed
    complex(r64), intent(in)    :: tx(nt), tnx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: sws(p*np), tpan(np+1)
    real(r64),    intent(inout) :: A(nt, np*p, M+1)
  end subroutine axissymlap_slpn_blockmat_nmode_r64

  subroutine axissymlap_dlp_blockmat_r64(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, A)
    integer(8),   intent(in)    :: nt, p, np, iside, iclosed
    complex(r64), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: sws(p*np), tpan(np+1)
    real(r64),    intent(inout) :: A(nt, np*p)
  end subroutine axissymlap_dlp_blockmat_r64

  subroutine axissymlap_dlp_blockmat_nmode_r64(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, A)
    integer(8),   intent(in)    :: nt, p, np, M, iside, iclosed
    complex(r64), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: sws(p*np), tpan(np+1)
    real(r64),    intent(inout) :: A(nt, np*p, M+1)
  end subroutine axissymlap_dlp_blockmat_nmode_r64

  subroutine axissymlap_dlpn_blockmat_r64(nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, A)
    integer(8),   intent(in)    :: nt, p, np, iside, iclosed
    complex(r64), intent(in)    :: tx(nt), tnx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: sws(p*np), tpan(np+1)
    real(r64),    intent(inout) :: A(nt, np*p)
  end subroutine axissymlap_dlpn_blockmat_r64

  subroutine axissymlap_dlpn_blockmat_nmode_r64(nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, A)
    integer(8),   intent(in)    :: nt, p, np, M, iside, iclosed
    complex(r64), intent(in)    :: tx(nt), tnx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: sws(p*np), tpan(np+1)
    real(r64),    intent(inout) :: A(nt, np*p, M+1)
  end subroutine axissymlap_dlpn_blockmat_nmode_r64

end module axissymlap_specialquad_mod
