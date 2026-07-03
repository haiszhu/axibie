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
