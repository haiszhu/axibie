! Top-level (free-standing) mex-facing wrappers for axissym_physop_sparse_mod.
! NOT inside a module -- module name mangling would break the mwrap binding.
! Prefix axps_ = AxiStokes + physop_sparse (axp_ is taken by physop).
!
! closesize / closeasm are SIGNATURES ONLY -- commented out until the module bodies land.
! Uncomment each wrapper together with its module procedure, then add the matching
! @function block to matlab/AxiStokes3D.mw.

subroutine axps_closesize_r64(nt, tx, t3dx, p, np, sx, sws, gate, tcxi, ntcx)
  use axissym_physop_sparse_mod, only: closesize => axissym_closesize_r64
  implicit none
  integer(8), intent(in)    :: nt, p, np
  complex(8), intent(in)    :: tx(nt), sx(p*np)
  real(8),    intent(in)    :: t3dx(3,nt), sws(p*np), gate
  real(8),    intent(inout) :: tcxi(np+1), ntcx
  call closesize(nt, tx, t3dx, p, np, sx, sws, gate, tcxi, ntcx)
end subroutine axps_closesize_r64

subroutine axps_closeslp_r64(nt, tx, t3dx, p, np, nang, sx, snx, sws, swxp, tpan, gate, &
                             s3dx, s3dnx, s3dw, pmodes, iside, iclosed, ntcx, tcxi, S_ij, idxall)
  use axissym_physop_sparse_mod, only: closeslp => axissym_closeslp_r64
  implicit none
  integer(8), intent(in)    :: nt, p, np, nang, pmodes, iside, iclosed, ntcx
  complex(8), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np)
  real(8),    intent(in)    :: t3dx(3,nt), sws(p*np), tpan(np+1), gate
  real(8),    intent(in)    :: s3dx(3, p*np*nang), s3dnx(3, p*np*nang), s3dw(p*np*nang), tcxi(np+1)
  real(8),    intent(inout) :: S_ij(ntcx, nang*p), idxall(ntcx)
  call closeslp(nt, tx, t3dx, p, np, nang, sx, snx, sws, swxp, tpan, gate, &
                s3dx, s3dnx, s3dw, pmodes, iside, iclosed, ntcx, tcxi, S_ij, idxall)
end subroutine axps_closeslp_r64

subroutine axps_closeslpn_r64(nt, tx, t3dx, tn3dx, p, np, nang, sx, snx, sws, swxp, tpan, gate, &
                              s3dx, s3dnx, s3dw, pmodes, iside, iclosed, ntcx, tcxi, S_ij, idxall)
  use axissym_physop_sparse_mod, only: closeslpn => axissym_closeslpn_r64
  implicit none
  integer(8), intent(in)    :: nt, p, np, nang, pmodes, iside, iclosed, ntcx
  complex(8), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np)
  real(8),    intent(in)    :: t3dx(3,nt), tn3dx(3,nt), sws(p*np), tpan(np+1), gate
  real(8),    intent(in)    :: s3dx(3, p*np*nang), s3dnx(3, p*np*nang), s3dw(p*np*nang), tcxi(np+1)
  real(8),    intent(inout) :: S_ij(ntcx, nang*p), idxall(ntcx)
  call closeslpn(nt, tx, t3dx, tn3dx, p, np, nang, sx, snx, sws, swxp, tpan, gate, &
                 s3dx, s3dnx, s3dw, pmodes, iside, iclosed, ntcx, tcxi, S_ij, idxall)
end subroutine axps_closeslpn_r64

subroutine axps_closedlp_r64(nt, tx, t3dx, p, np, nang, sx, snx, sws, swxp, tpan, gate, &
                             s3dx, s3dnx, s3dw, pmodes, iside, iclosed, ntcx, tcxi, S_ij, idxall)
  use axissym_physop_sparse_mod, only: closedlp => axissym_closedlp_r64
  implicit none
  integer(8), intent(in)    :: nt, p, np, nang, pmodes, iside, iclosed, ntcx
  complex(8), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np)
  real(8),    intent(in)    :: t3dx(3,nt), sws(p*np), tpan(np+1), gate
  real(8),    intent(in)    :: s3dx(3, p*np*nang), s3dnx(3, p*np*nang), s3dw(p*np*nang), tcxi(np+1)
  real(8),    intent(inout) :: S_ij(ntcx, nang*p), idxall(ntcx)
  call closedlp(nt, tx, t3dx, p, np, nang, sx, snx, sws, swxp, tpan, gate, &
                s3dx, s3dnx, s3dw, pmodes, iside, iclosed, ntcx, tcxi, S_ij, idxall)
end subroutine axps_closedlp_r64

subroutine axps_closeasm_r64(iphys, ilayer, nc, nt, tx, t3dx, p, np, nang, s3dx, s3dnx, s3dw, &
                             mu, ntcx, tcxi, idxall, S_ij, iself, A)
  use axissym_physop_sparse_mod, only: closeasm => axissym_closeasm_r64
  implicit none
  integer(8), intent(in)    :: iphys, ilayer, nc, nt, p, np, nang, ntcx, iself
  complex(8), intent(in)    :: tx(nt)
  real(8),    intent(in)    :: t3dx(3,nt), s3dx(3,p*np*nang), s3dnx(3,p*np*nang), s3dw(p*np*nang), mu
  real(8),    intent(in)    :: tcxi(np+1), idxall(ntcx), S_ij(ntcx, nang*p)
  real(8),    intent(inout) :: A(*)
  call closeasm(iphys, ilayer, nc, nt, tx, t3dx, p, np, nang, s3dx, s3dnx, s3dw, &
                mu, ntcx, tcxi, idxall, S_ij, iself, A)
end subroutine axps_closeasm_r64

subroutine axps_closecorrapply_r64(nc, nt, p, np, nang, ntcx, tcxi, idxall, S_ij, sig, iself, u)
  use axissym_physop_sparse_mod, only: corrapply => axissym_closecorrapply_r64
  implicit none
  integer(8), intent(in)    :: nc, nt, p, np, nang, ntcx, iself
  real(8),    intent(in)    :: tcxi(np+1), idxall(ntcx), S_ij(nc*ntcx, nc*nang*p), sig(nc*p*np*nang)
  real(8),    intent(inout) :: u(*)
  call corrapply(nc, nt, p, np, nang, ntcx, tcxi, idxall, S_ij, sig, iself, u)
end subroutine axps_closecorrapply_r64

subroutine axps_closestokslp_r64(nt, tx, t3dx, p, np, nang, sx, snx, sws, swxp, tpan, gate, &
                                 s3dx, s3dnx, s3dw, pmodes, iside, iclosed, mu, ntcx, tcxi, S_ij, idxall)
  use axissym_physop_sparse_mod, only: closestokslp => axissym_closestokslp_r64
  implicit none
  integer(8), intent(in)    :: nt, p, np, nang, pmodes, iside, iclosed, ntcx
  complex(8), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np)
  real(8),    intent(in)    :: t3dx(3,nt), sws(p*np), tpan(np+1), gate, mu
  real(8),    intent(in)    :: s3dx(3, p*np*nang), s3dnx(3, p*np*nang), s3dw(p*np*nang), tcxi(np+1)
  real(8),    intent(inout) :: S_ij(3*ntcx, 3*nang*p), idxall(ntcx)
  call closestokslp(nt, tx, t3dx, p, np, nang, sx, snx, sws, swxp, tpan, gate, &
                    s3dx, s3dnx, s3dw, pmodes, iside, iclosed, mu, ntcx, tcxi, S_ij, idxall)
end subroutine axps_closestokslp_r64

subroutine axps_closestokslpn_r64(nt, tx, t3dx, tn3dx, p, np, nang, sx, snx, sws, swxp, tpan, gate, &
                                  s3dx, s3dnx, s3dw, pmodes, iside, iclosed, mu, ntcx, tcxi, S_ij, idxall)
  use axissym_physop_sparse_mod, only: closestokslpn => axissym_closestokslpn_r64
  implicit none
  integer(8), intent(in)    :: nt, p, np, nang, pmodes, iside, iclosed, ntcx
  complex(8), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np)
  real(8),    intent(in)    :: t3dx(3,nt), tn3dx(3,nt), sws(p*np), tpan(np+1), gate, mu
  real(8),    intent(in)    :: s3dx(3, p*np*nang), s3dnx(3, p*np*nang), s3dw(p*np*nang), tcxi(np+1)
  real(8),    intent(inout) :: S_ij(3*ntcx, 3*nang*p), idxall(ntcx)
  call closestokslpn(nt, tx, t3dx, tn3dx, p, np, nang, sx, snx, sws, swxp, tpan, gate, &
                     s3dx, s3dnx, s3dw, pmodes, iside, iclosed, mu, ntcx, tcxi, S_ij, idxall)
end subroutine axps_closestokslpn_r64

subroutine axps_closestokdlpn_r64(nt, tx, t3dx, tn3dx, p, np, nang, sx, snx, sws, swxp, tpan, gate, &
                                  s3dx, s3dnx, s3dw, pmodes, iside, iclosed, mu, ntcx, tcxi, S_ij, idxall)
  use axissym_physop_sparse_mod, only: closestokdlpn => axissym_closestokdlpn_r64
  implicit none
  integer(8), intent(in)    :: nt, p, np, nang, pmodes, iside, iclosed, ntcx
  complex(8), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np)
  real(8),    intent(in)    :: t3dx(3,nt), tn3dx(3,nt), sws(p*np), tpan(np+1), gate, mu
  real(8),    intent(in)    :: s3dx(3, p*np*nang), s3dnx(3, p*np*nang), s3dw(p*np*nang), tcxi(np+1)
  real(8),    intent(inout) :: S_ij(3*ntcx, 3*nang*p), idxall(ntcx)
  call closestokdlpn(nt, tx, t3dx, tn3dx, p, np, nang, sx, snx, sws, swxp, tpan, gate, &
                     s3dx, s3dnx, s3dw, pmodes, iside, iclosed, mu, ntcx, tcxi, S_ij, idxall)
end subroutine axps_closestokdlpn_r64

subroutine axps_closelapsdlp_panel_r64(nt, tx, t3dx, tn3dx, p, nang, sx, snx, sws, swxp, sxlo, sxhi, tpan, gate, &
                             s3dx, s3dnx, s3dw, pmodes, iside, iclosed, ntcx, tcxi, idxall, iform, As3d, Ad3d)
  use axissym_physop_sparse_mod, only: panel => axissym_closelapsdlp_panel_r64
  implicit none
  integer(8), intent(in)    :: nt, p, nang, pmodes, iside, iclosed, ntcx, iform
  complex(8), intent(in)    :: tx(nt), sx(p), snx(p), swxp(p), sxlo, sxhi
  real(8),    intent(in)    :: t3dx(3,nt), tn3dx(3,nt), sws(p), tpan(*), gate, s3dx(3,*), s3dnx(3,*), s3dw(*), tcxi(*), idxall(*)
  real(8),    intent(inout) :: As3d(nt,nang*p), Ad3d(nt,nang*p)
  call panel(nt, tx, t3dx, tn3dx, p, nang, sx, snx, sws, swxp, sxlo, sxhi, tpan, gate, &
             s3dx, s3dnx, s3dw, pmodes, iside, iclosed, ntcx, tcxi, idxall, iform, As3d, Ad3d)
end subroutine axps_closelapsdlp_panel_r64

subroutine axps_closestoksdlp_panel_r64(nt, tx, t3dx, tn3dx, p, nang, sx, snx, sws, swxp, sxlo, sxhi, tpan, gate, &
                             s3dx, s3dnx, s3dw, pmodes, iside, iclosed, mu, ntcx, tcxi, idxall, iform, As3d, Ad3d)
  use axissym_physop_sparse_mod, only: panel => axissym_closestoksdlp_panel_r64
  implicit none
  integer(8), intent(in)    :: nt, p, nang, pmodes, iside, iclosed, ntcx, iform
  complex(8), intent(in)    :: tx(nt), sx(p), snx(p), swxp(p), sxlo, sxhi
  real(8),    intent(in)    :: t3dx(3,nt), tn3dx(3,nt), sws(p), tpan(*), gate, s3dx(3,*), s3dnx(3,*), s3dw(*), mu, tcxi(*), idxall(*)
  real(8),    intent(inout) :: As3d(3*nt,3*nang*p), Ad3d(3*nt,3*nang*p)
  call panel(nt, tx, t3dx, tn3dx, p, nang, sx, snx, sws, swxp, sxlo, sxhi, tpan, gate, &
             s3dx, s3dnx, s3dw, pmodes, iside, iclosed, mu, ntcx, tcxi, idxall, iform, As3d, Ad3d)
end subroutine axps_closestoksdlp_panel_r64

subroutine axps_naivelapsdlp_physmat_r64(nt, tx, t3dx, p, np, nang, sx, snx, sws, pmodes, ifself, nrA, As3d, Ad3d)
  use axissym_physop_sparse_mod, only: physmat => axissym_naivelapsdlp_physmat_r64
  implicit none
  integer(8), intent(in)    :: nt, p, np, nang, pmodes, ifself, nrA
  complex(8), intent(in)    :: tx(nt), sx(p*np), snx(p*np)
  real(8),    intent(in)    :: t3dx(3,nt), sws(p*np)
  real(8),    intent(inout) :: As3d(nrA,np*p*nang), Ad3d(nrA,np*p*nang)
  call physmat(nt, tx, t3dx, p, np, nang, sx, snx, sws, pmodes, ifself, nrA, As3d, Ad3d)
end subroutine axps_naivelapsdlp_physmat_r64
