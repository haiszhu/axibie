! Top-level (free-standing) mex-facing wrappers for axissym_physop_mod.
! NOT inside a module -- module name mangling would break the mwrap binding.

subroutine axp_physmat_r64(iphys, ilayer, nc, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, &
                           M, iside, iclosed, mu, R, &
                           iAb, jAb, xAb, iAi, jAi, xAi, iF, jF, xF, iFi, jFi, xFi, iT, jT, xT)
  use axissym_physop_mod, only: physmat => axissym_physmat_r64
  implicit none
  integer(8), intent(in)    :: iphys, ilayer, nc, p, np, M, iside, iclosed
  complex(8), intent(in)    :: sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
  real(8),    intent(in)    :: sws(p*np), tpan(np+1), mu, R(3,3)
  real(8),    intent(inout) :: iAb(nc*p*np*(2*M+1)+1), jAb((2*M+1)*(nc*p*np)**2)
  complex(8), intent(inout) :: xAb((2*M+1)*(nc*p*np)**2)
  real(8),    intent(inout) :: iAi(nc*p*np*(2*M+1)+1), jAi((2*M+1)*(nc*p*np)**2)
  complex(8), intent(inout) :: xAi((2*M+1)*(nc*p*np)**2)
  real(8),    intent(inout) :: iF(nc*p*np*(2*M+1)+1), jF((2*M+1)**2*(nc*p*np))
  complex(8), intent(inout) :: xF((2*M+1)**2*(nc*p*np))
  real(8),    intent(inout) :: iFi(nc*p*np*(2*M+1)+1), jFi((2*M+1)**2*(nc*p*np))
  complex(8), intent(inout) :: xFi((2*M+1)**2*(nc*p*np))
  real(8),    intent(inout) :: iT(nc*p*np*(2*M+1)+1), jT(nc*nc*p*np*(2*M+1))
  real(8),    intent(inout) :: xT(nc*nc*p*np*(2*M+1))
  call physmat(iphys, ilayer, nc, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, &
               M, iside, iclosed, mu, R, &
               iAb, jAb, xAb, iAi, jAi, xAi, iF, jF, xF, iFi, jFi, xFi, iT, jT, xT)
end subroutine axp_physmat_r64

subroutine axp_offdiagphysmat_r64(iphys, ilayer, nc, nt, tx, tphi, tnx, tnphi, Ct, Rt, &
                                  p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, Cs, Rs, &
                                  M, iside, iclosed, mu, iforce, &
                                  B, iF, jF, xF, near)
  use axissym_physop_mod, only: offdiag => axissym_offdiagphysmat_r64
  implicit none
  integer(8), intent(in)    :: iphys, ilayer, nc, nt, p, np, M, iside, iclosed, iforce
  complex(8), intent(in)    :: tx(nt), tnx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
  real(8),    intent(in)    :: tphi(nt), tnphi(nt), Ct(3), Rt(3,3), Cs(3), Rs(3,3), sws(p*np), tpan(np+1), mu
  complex(8), intent(inout) :: B(nc*nt, nc*p*np*(2*M+1))
  real(8),    intent(inout) :: iF(nc*p*np*(2*M+1)+1), jF((2*M+1)**2*(nc*p*np))
  complex(8), intent(inout) :: xF((2*M+1)**2*(nc*p*np))
  real(8),    intent(inout) :: near(nt)
  call offdiag(iphys, ilayer, nc, nt, tx, tphi, tnx, tnphi, Ct, Rt, &
               p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, Cs, Rs, &
               M, iside, iclosed, mu, iforce, B, iF, jF, xF, near)
end subroutine axp_offdiagphysmat_r64
