! Top-level (free-standing) mex-facing wrappers for axissym_physop_mod.
! NOT inside a module -- module name mangling would break the mwrap binding.

subroutine axp_physmat_r64(iphys, ilayer, nc, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, &
                           M, iside, iclosed, mu, R, xA, &
                           iAb, jAb, xAb, iAi, jAi, xAi, iF, jF, xF, iFi, jFi, xFi, iT, jT, xT)
  use axissym_physop_mod, only: physmat => axissym_physmat_r64
  implicit none
  integer(8), intent(in)    :: iphys, ilayer, nc, p, np, M, iside, iclosed
  complex(8), intent(in)    :: sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
  real(8),    intent(in)    :: sws(p*np), tpan(np+1), mu, R(3,3)
  real(8),    intent(inout) :: xA(nc*p*np*(2*M+1), nc*p*np*(2*M+1))
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
               M, iside, iclosed, mu, R, xA, &
               iAb, jAb, xAb, iAi, jAi, xAi, iF, jF, xF, iFi, jFi, xFi, iT, jT, xT)
end subroutine axp_physmat_r64

subroutine axp_offdiagphysmat_r64(iphys, ilayer, nc, nt, ns, Xt, Ntg, Ct, Rt, &
                                  p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, Cs, Rs, &
                                  M, iside, iclosed, mu, iforce, &
                                  At, B, iF, jF, xF, near)
  use axissym_physop_mod, only: offdiag => axissym_offdiagphysmat_r64
  implicit none
  integer(8), intent(in)    :: iphys, ilayer, nc, nt, ns, p, np, M, iside, iclosed, iforce
  complex(8), intent(in)    :: sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
  real(8),    intent(in)    :: Xt(3,nt), Ntg(3,nt), Ct(3), Rt(3,3), Cs(3), Rs(3,3), sws(p*np), tpan(np+1), mu
  real(8),    intent(inout) :: At(nc*nt, nc*ns)
  complex(8), intent(inout) :: B(nc*nt, nc*p*np*(2*M+1))
  real(8),    intent(inout) :: iF(nc*p*np*(2*M+1)+1), jF((2*M+1)**2*(nc*p*np))
  complex(8), intent(inout) :: xF((2*M+1)**2*(nc*p*np))
  real(8),    intent(inout) :: near(nt)
  call offdiag(iphys, ilayer, nc, nt, ns, Xt, Ntg, Ct, Rt, &
               p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, Cs, Rs, &
               M, iside, iclosed, mu, iforce, At, B, iF, jF, xF, near)
end subroutine axp_offdiagphysmat_r64

subroutine axp_physmat_setup_r64(ikernel, ilayer, params, iinter, K, p, np, pmodes, iside, iclosed, &
                                 geomoff, tpanoff, nsx, ntpan, sx, snx, sws, swxp, tpan, &
                                 ntarg, targoff, targ, targnx, Rm, Cc, nrA, ncA, A)
  use axissym_physop_mod, only: physmatsetup => axissym_physmat_setup
  implicit none
  integer(8), intent(in)    :: ikernel, ilayer, iinter, K, iside, iclosed, nsx, ntpan, ntarg, nrA, ncA
  integer(8), intent(in)    :: p(K), np(K), pmodes(K), geomoff(K+1), tpanoff(K+1), targoff(K+1)
  complex(8), intent(in)    :: params, sx(nsx), snx(nsx), swxp(nsx)
  real(8),    intent(in)    :: sws(nsx), tpan(ntpan), targ(3,ntarg), targnx(3,ntarg), Rm(3,3,K), Cc(3,K)
  real(8),    intent(inout) :: A(nrA, ncA)                       ! caller pre-zeroed (untouched blocks stay 0)
  call physmatsetup(ikernel, ilayer, params, iinter, K, p, np, pmodes, iside, iclosed, &
                    geomoff, tpanoff, nsx, ntpan, sx, snx, sws, swxp, tpan, &
                    ntarg, targoff, targ, targnx, Rm, Cc, nrA, ncA, A)
end subroutine axp_physmat_setup_r64

subroutine axp_modemat_setup_r64(ikernel, ilayer, params, iinter, K, p, np, pmodes, iside, iclosed, &
                                 geomoff, tpanoff, nsx, ntpan, sx, snx, sws, swxp, tpan, &
                                 ntarg, targoff, targ, targnx, Rm, Cc, nrA, ncA, A)
  use axissym_physop_mod, only: modemat => axissym_modemat_setup
  implicit none
  integer(8), intent(in)    :: ikernel, ilayer, iinter, K, iside, iclosed, nsx, ntpan, ntarg, nrA, ncA
  integer(8), intent(in)    :: p(K), np(K), pmodes(K), geomoff(K+1), tpanoff(K+1), targoff(K+1)
  complex(8), intent(in)    :: params, sx(nsx), snx(nsx), swxp(nsx)
  real(8),    intent(in)    :: sws(nsx), tpan(ntpan), targ(3,ntarg), targnx(3,ntarg), Rm(3,3,K), Cc(3,K)
  complex(8), intent(inout) :: A(nrA, ncA, *)                    ! 3rd dim = maxval(pmodes)+1, caller pre-zeroed
  call modemat(ikernel, ilayer, params, iinter, K, p, np, pmodes, iside, iclosed, &
               geomoff, tpanoff, nsx, ntpan, sx, snx, sws, swxp, tpan, &
               ntarg, targoff, targ, targnx, Rm, Cc, nrA, ncA, A)
end subroutine axp_modemat_setup_r64
