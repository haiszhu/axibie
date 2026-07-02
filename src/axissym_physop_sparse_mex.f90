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
                             ntcx, tcxi, idxall, S_ij, iself, A)
  use axissym_physop_sparse_mod, only: closeasm => axissym_closeasm_r64
  implicit none
  integer(8), intent(in)    :: iphys, ilayer, nc, nt, p, np, nang, ntcx, iself
  complex(8), intent(in)    :: tx(nt)
  real(8),    intent(in)    :: t3dx(3,nt), s3dx(3,p*np*nang), s3dnx(3,p*np*nang), s3dw(p*np*nang)
  real(8),    intent(in)    :: tcxi(np+1), idxall(ntcx), S_ij(ntcx, nang*p)
  real(8),    intent(inout) :: A(*)
  call closeasm(iphys, ilayer, nc, nt, tx, t3dx, p, np, nang, s3dx, s3dnx, s3dw, &
                ntcx, tcxi, idxall, S_ij, iself, A)
end subroutine axps_closeasm_r64

subroutine axps_closecorrapply_r64(nt, p, np, nang, ntcx, tcxi, idxall, S_ij, sig, iself, u)
  use axissym_physop_sparse_mod, only: corrapply => axissym_closecorrapply_r64
  implicit none
  integer(8), intent(in)    :: nt, p, np, nang, ntcx, iself
  real(8),    intent(in)    :: tcxi(np+1), idxall(ntcx), S_ij(ntcx, nang*p), sig(p*np*nang)
  real(8),    intent(inout) :: u(*)
  call corrapply(nt, p, np, nang, ntcx, tcxi, idxall, S_ij, sig, iself, u)
end subroutine axps_closecorrapply_r64
