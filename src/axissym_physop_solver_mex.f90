! Top-level (free-standing) mex-facing wrappers for axissym_physop_solver_mod.
! NOT inside a module -- module name mangling would break the mwrap binding.
! Prefix axpso_ = AxiStokes + physop_solver (axps_ = physop_sparse is taken;
! collision rule extends the module letter to two).

subroutine axpso_stokslpn_release_r64(ok)
  use axissym_physop_solver_mod, only: release => axissym_stokslpn_release_r64
  implicit none
  real(8), intent(inout) :: ok
  call release(ok)
end subroutine axpso_stokslpn_release_r64

subroutine axpso_stokslpn_selfsetup_r64(K, p, np, nang, pmodes, iside, iclosed, gate, mu, &
                                        sx, snx, sws, swxp, tpan, ntcxs, nbytes)
  use axissym_physop_solver_mod, only: selfsetup => axissym_stokslpn_selfsetup_r64
  implicit none
  integer(8), intent(in)    :: K, p, np, nang, pmodes, iside, iclosed
  real(8),    intent(in)    :: gate, mu, sws(p*np,K), tpan(np+1,K)
  complex(8), intent(in)    :: sx(p*np,K), snx(p*np,K), swxp(p*np,K)
  real(8),    intent(inout) :: ntcxs(K), nbytes
  call selfsetup(K, p, np, nang, pmodes, iside, iclosed, gate, mu, sx, snx, sws, swxp, tpan, ntcxs, nbytes)
end subroutine axpso_stokslpn_selfsetup_r64

subroutine axpso_stokslpn_getself_r64(k, ntcxk, p, np, nang, Sk, ik, tcxik)
  use axissym_physop_solver_mod, only: getself => axissym_stokslpn_getself_r64
  implicit none
  integer(8), intent(in)    :: k, ntcxk, p, np, nang
  real(8),    intent(inout) :: Sk(3*ntcxk, 3*nang*p), ik(ntcxk), tcxik(np+1)
  call getself(k, ntcxk, p, np, nang, Sk, ik, tcxik)
end subroutine axpso_stokslpn_getself_r64

subroutine axpso_stokslpn_crosssetup_r64(K, p, np, nang, pmodes, iside, iclosed, gate, mu, rball, &
                                         Xall, Nall, Rm, Cc, sx, snx, sws, swxp, tpan, &
                                         ntcxx, nuo, npair, nbytes)
  use axissym_physop_solver_mod, only: crosssetup => axissym_stokslpn_crosssetup_r64
  implicit none
  integer(8), intent(in)    :: K, p, np, nang, pmodes, iside, iclosed
  real(8),    intent(in)    :: gate, mu, rball, Xall(3,K*p*np*nang), Nall(3,K*p*np*nang)
  real(8),    intent(in)    :: Rm(3,3,K), Cc(3,K), sws(p*np,K), tpan(np+1,K)
  complex(8), intent(in)    :: sx(p*np,K), snx(p*np,K), swxp(p*np,K)
  real(8),    intent(inout) :: ntcxx(K), nuo(K), npair, nbytes
  call crosssetup(K, p, np, nang, pmodes, iside, iclosed, gate, mu, rball, &
                  Xall, Nall, Rm, Cc, sx, snx, sws, swxp, tpan, ntcxx, nuo, npair, nbytes)
end subroutine axpso_stokslpn_crosssetup_r64

subroutine axpso_stokslpn_getcross_r64(k, ntcxk, nuk, p, np, nang, Sk, idxk, tcxik, canonk)
  use axissym_physop_solver_mod, only: getcross => axissym_stokslpn_getcross_r64
  implicit none
  integer(8), intent(in)    :: k, ntcxk, nuk, p, np, nang
  real(8),    intent(inout) :: Sk(3*ntcxk, 3*nang*p), idxk(ntcxk), tcxik(np+1), canonk(nuk)
  call getcross(k, ntcxk, nuk, p, np, nang, Sk, idxk, tcxik, canonk)
end subroutine axpso_stokslpn_getcross_r64

subroutine axpso_stokslpn_corrapply_r64(K, n, Rm, x, u)
  use axissym_physop_solver_mod, only: corrapply => axissym_stokslpn_corrapply_r64
  implicit none
  integer(8), intent(in)    :: K, n
  real(8),    intent(in)    :: Rm(3,3,K), x(n)
  real(8),    intent(inout) :: u(n)
  call corrapply(K, n, Rm, x, u)
end subroutine axpso_stokslpn_corrapply_r64

subroutine axpso_stokslp_fieldeval_r64(K, Me, p, np, nang, pmodes, iside, iclosed, gate, mu, rball, &
                                       Pe, sigma, Rm, Cc, sx, snx, sws, swxp, tpan, u)
  use axissym_physop_solver_mod, only: fieldeval => axissym_stokslp_fieldeval_r64
  implicit none
  integer(8), intent(in)    :: K, Me, p, np, nang, pmodes, iside, iclosed
  real(8),    intent(in)    :: gate, mu, rball, Pe(3,Me), sigma(3*K*p*np*nang)
  real(8),    intent(in)    :: Rm(3,3,K), Cc(3,K), sws(p*np,K), tpan(np+1,K)
  complex(8), intent(in)    :: sx(p*np,K), snx(p*np,K), swxp(p*np,K)
  real(8),    intent(inout) :: u(3*Me)
  call fieldeval(K, Me, p, np, nang, pmodes, iside, iclosed, gate, mu, rball, &
                 Pe, sigma, Rm, Cc, sx, snx, sws, swxp, tpan, u)
end subroutine axpso_stokslp_fieldeval_r64

subroutine axpso_lapslpn_release_r64(ok)
  use axissym_physop_solver_mod, only: release => axissym_lapslpn_release_r64
  implicit none
  real(8), intent(inout) :: ok
  call release(ok)
end subroutine axpso_lapslpn_release_r64

subroutine axpso_lapslpn_selfsetup_r64(K, p, np, nang, pmodes, iside, iclosed, gate, &
                                       sx, snx, sws, swxp, tpan, ntcxs, nbytes)
  use axissym_physop_solver_mod, only: selfsetup => axissym_lapslpn_selfsetup_r64
  implicit none
  integer(8), intent(in)    :: K, p, np, nang, pmodes, iside, iclosed
  real(8),    intent(in)    :: gate, sws(p*np,K), tpan(np+1,K)
  complex(8), intent(in)    :: sx(p*np,K), snx(p*np,K), swxp(p*np,K)
  real(8),    intent(inout) :: ntcxs(K), nbytes
  call selfsetup(K, p, np, nang, pmodes, iside, iclosed, gate, sx, snx, sws, swxp, tpan, ntcxs, nbytes)
end subroutine axpso_lapslpn_selfsetup_r64

subroutine axpso_lapslpn_getself_r64(k, ntcxk, p, np, nang, Sk, ik, tcxik)
  use axissym_physop_solver_mod, only: getself => axissym_lapslpn_getself_r64
  implicit none
  integer(8), intent(in)    :: k, ntcxk, p, np, nang
  real(8),    intent(inout) :: Sk(ntcxk, nang*p), ik(ntcxk), tcxik(np+1)
  call getself(k, ntcxk, p, np, nang, Sk, ik, tcxik)
end subroutine axpso_lapslpn_getself_r64

subroutine axpso_lapslpn_crosssetup_r64(K, p, np, nang, pmodes, iside, iclosed, gate, rball, &
                                        Xall, Nall, Rm, Cc, sx, snx, sws, swxp, tpan, &
                                        ntcxx, nuo, npair, nbytes)
  use axissym_physop_solver_mod, only: crosssetup => axissym_lapslpn_crosssetup_r64
  implicit none
  integer(8), intent(in)    :: K, p, np, nang, pmodes, iside, iclosed
  real(8),    intent(in)    :: gate, rball, Xall(3,K*p*np*nang), Nall(3,K*p*np*nang)
  real(8),    intent(in)    :: Rm(3,3,K), Cc(3,K), sws(p*np,K), tpan(np+1,K)
  complex(8), intent(in)    :: sx(p*np,K), snx(p*np,K), swxp(p*np,K)
  real(8),    intent(inout) :: ntcxx(K), nuo(K), npair, nbytes
  call crosssetup(K, p, np, nang, pmodes, iside, iclosed, gate, rball, &
                  Xall, Nall, Rm, Cc, sx, snx, sws, swxp, tpan, ntcxx, nuo, npair, nbytes)
end subroutine axpso_lapslpn_crosssetup_r64

subroutine axpso_lapslpn_getcross_r64(k, ntcxk, nuk, p, np, nang, Sk, idxk, tcxik, canonk)
  use axissym_physop_solver_mod, only: getcross => axissym_lapslpn_getcross_r64
  implicit none
  integer(8), intent(in)    :: k, ntcxk, nuk, p, np, nang
  real(8),    intent(inout) :: Sk(ntcxk, nang*p), idxk(ntcxk), tcxik(np+1), canonk(nuk)
  call getcross(k, ntcxk, nuk, p, np, nang, Sk, idxk, tcxik, canonk)
end subroutine axpso_lapslpn_getcross_r64

subroutine axpso_lapslpn_corrapply_r64(K, n, x, u)
  use axissym_physop_solver_mod, only: corrapply => axissym_lapslpn_corrapply_r64
  implicit none
  integer(8), intent(in)    :: K, n
  real(8),    intent(in)    :: x(n)
  real(8),    intent(inout) :: u(n)
  call corrapply(K, n, x, u)
end subroutine axpso_lapslpn_corrapply_r64

! ==== WIP: handle-based near-correction API (BOUNDARY, mwrap-facing, intrinsics only) ====
! Builds the descriptor from scalars, runs corr_setup, stores the resulting sparse_corr_t in
! the module registry, returns its integer handle (as a double). (empty skeleton -- for review)
subroutine axpso_corr_setup_r64(ikernel, ilayer, ipanel, dbg, params, iinter, K, p, np, pmodes, iside, iclosed, gate, &
                                geomoff, tpanoff, nsx, ntpan, sx, snx, sws, swxp, tpan, &
                                rball, ntarg, targoff, targ, targnx, Rm, Cc, &
                                handle, nbytes)
  use axissym_physop_solver_mod, only: axissym_corr_setup
  implicit none
  integer(8), intent(in)    :: ikernel, ilayer, ipanel, dbg, iinter, K, iside, iclosed, ntarg, nsx, ntpan
  integer(8), intent(in)    :: p(K), np(K), pmodes(K), geomoff(K+1), tpanoff(K+1), targoff(K+1)
  complex(8), intent(in)    :: params
  real(8),    intent(in)    :: gate, rball
  complex(8), intent(in)    :: sx(nsx), snx(nsx), swxp(nsx)
  real(8),    intent(in)    :: sws(nsx), tpan(ntpan)
  real(8),    intent(in)    :: targ(3,ntarg), targnx(3,ntarg), Rm(3,3,K), Cc(3,K)
  real(8),    intent(inout) :: handle, nbytes
  integer(8) :: h, nb
  h = nint(handle, 8)                                ! 0 = monolithic ; > 0 = handle from axpso_create_mex
  call axissym_corr_setup(ikernel, ilayer, ipanel, dbg, params, iinter, K, p, np, pmodes, iside, iclosed, gate, &
                          geomoff, tpanoff, nsx, ntpan, sx, snx, sws, swxp, tpan, &
                          rball, ntarg, targoff, targ, targnx, Rm, Cc, h, nb)
  handle = real(h, 8); nbytes = real(nb, 8)
end subroutine axpso_corr_setup_r64

subroutine axpso_close_setup_r64(ikernel, ilayer, ipanel, dbg, params, iinter, K, p, np, pmodes, iside, iclosed, gate, &
                                geomoff, tpanoff, nsx, ntpan, sx, snx, sws, swxp, tpan, &
                                rball, ntarg, targoff, targ, targnx, Rm, Cc, &
                                handle, nbytes)
  use axissym_physop_solver_mod, only: axissym_close_setup
  implicit none
  integer(8), intent(in)    :: ikernel, ilayer, ipanel, dbg, iinter, K, iside, iclosed, ntarg, nsx, ntpan
  integer(8), intent(in)    :: p(K), np(K), pmodes(K), geomoff(K+1), tpanoff(K+1), targoff(K+1)
  complex(8), intent(in)    :: params
  real(8),    intent(in)    :: gate, rball
  complex(8), intent(in)    :: sx(nsx), snx(nsx), swxp(nsx)
  real(8),    intent(in)    :: sws(nsx), tpan(ntpan)
  real(8),    intent(in)    :: targ(3,ntarg), targnx(3,ntarg), Rm(3,3,K), Cc(3,K)
  real(8),    intent(inout) :: handle, nbytes
  integer(8) :: h, nb
  h = nint(handle, 8)                                ! 0 = monolithic ; > 0 = handle from axpso_create_mex
  call axissym_close_setup(ikernel, ilayer, ipanel, dbg, params, iinter, K, p, np, pmodes, iside, iclosed, gate, &
                          geomoff, tpanoff, nsx, ntpan, sx, snx, sws, swxp, tpan, &
                          rball, ntarg, targoff, targ, targnx, Rm, Cc, h, nb)
  handle = real(h, 8); nbytes = real(nb, 8)
end subroutine axpso_close_setup_r64

subroutine axpso_create_r64(ikernel, ilayer, params, iinter, K, p, np, pmodes, iside, iclosed, gate, &
                                geomoff, tpanoff, nsx, ntpan, sx, snx, sws, swxp, tpan, &
                                rball, ntarg, targoff, targ, targnx, Rm, Cc, &
                                handle, nbytes)
  use axissym_physop_solver_mod, only: axissym_create
  implicit none
  integer(8), intent(in)    :: ikernel, ilayer, iinter, K, iside, iclosed, ntarg, nsx, ntpan
  integer(8), intent(in)    :: p(K), np(K), pmodes(K), geomoff(K+1), tpanoff(K+1), targoff(K+1)
  complex(8), intent(in)    :: params
  real(8),    intent(in)    :: gate, rball
  complex(8), intent(in)    :: sx(nsx), snx(nsx), swxp(nsx)
  real(8),    intent(in)    :: sws(nsx), tpan(ntpan)
  real(8),    intent(in)    :: targ(3,ntarg), targnx(3,ntarg), Rm(3,3,K), Cc(3,K)
  real(8),    intent(inout) :: handle, nbytes
  integer(8) :: h, nb
  call axissym_create(ikernel, ilayer, params, iinter, K, p, np, pmodes, iside, iclosed, gate, &
                          geomoff, tpanoff, nsx, ntpan, sx, snx, sws, swxp, tpan, &
                          rball, ntarg, targoff, targ, targnx, Rm, Cc, h, nb)
  handle = real(h, 8); nbytes = real(nb, 8)
end subroutine axpso_create_r64

subroutine axpso_corr_get_r64(handle, k, ntcxk, nc, p, np, pmodes, blk, ik, tcxik)
  use axissym_physop_solver_mod, only: axissym_corr_get
  implicit none
  real(8),    intent(in)    :: handle
  integer(8), intent(in)    :: k, ntcxk, nc, p, np, pmodes
  real(8),    intent(inout) :: blk(nc*nc*ntcxk*(2*pmodes+1)*p), ik(ntcxk), tcxik(np+1)
  call axissym_corr_get(nint(handle,8), k, ntcxk, nc, p, np, pmodes, blk, ik, tcxik)
end subroutine axpso_corr_get_r64

subroutine axpso_corr_set_r64(handle, nvals, vals, niks, iks)
  use axissym_physop_solver_mod, only: axissym_corr_set
  implicit none
  real(8),    intent(in) :: handle
  integer(8), intent(in) :: nvals, niks
  real(8),    intent(in) :: vals(nvals), iks(niks)
  call axissym_corr_set(nint(handle,8), nvals, vals, niks, iks)
end subroutine axpso_corr_set_r64

subroutine axpso_corr2dense_get_r64(handle, nrA, ncA, A)
  use axissym_physop_solver_mod, only: axissym_corr2dense_get
  implicit none
  real(8),    intent(in)    :: handle
  integer(8), intent(in)    :: nrA, ncA
  real(8),    intent(inout) :: A(nrA, ncA)
  call axissym_corr2dense_get(nint(handle,8), nrA, ncA, A)
end subroutine axpso_corr2dense_get_r64

subroutine axpso_close2corr_set_r64(handle, K, geomoff, nsx, sx, snx, sws, targoff, ntarg, targ, targnx, Cc)
  use axissym_physop_solver_mod, only: axissym_close2corr_set
  implicit none
  real(8),    intent(in) :: handle
  integer(8), intent(in) :: K, nsx, ntarg, geomoff(K+1), targoff(K+1)
  complex(8), intent(in) :: sx(nsx), snx(nsx)
  real(8),    intent(in) :: sws(nsx), targ(3,ntarg), targnx(3,ntarg), Cc(3,K)
  call axissym_close2corr_set(nint(handle,8), K, geomoff, nsx, sx, snx, sws, targoff, ntarg, targ, targnx, Cc)
end subroutine axpso_close2corr_set_r64

subroutine axpso_corr_apply_r64(handle, nx, x, nu, u)
  use axissym_physop_solver_mod, only: axissym_corr_apply
  implicit none
  real(8),    intent(in)    :: handle
  integer(8), intent(in)    :: nx, nu
  real(8),    intent(in)    :: x(nx)
  real(8),    intent(inout) :: u(nu)
  call axissym_corr_apply(nint(handle,8), nx, x, nu, u)
end subroutine axpso_corr_apply_r64

subroutine axpso_corr_free_r64(handle)
  use axissym_physop_solver_mod, only: axissym_corr_free
  implicit none
  real(8), intent(in) :: handle
  call axissym_corr_free(nint(handle,8))
end subroutine axpso_corr_free_r64
