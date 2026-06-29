! kdtree_mex.f90 — top-level mwrap wrappers (axk_ = ax + kdtree) over kdtree_mod.
! Stateful build / ball / free sequence; backend (KD_NBODY=1, KD_TAIYA=2) is required.

subroutine axk_kdtree_build_r64(backend, np, pts, leafsize)
  use kdtree_mod, only: kdtree_build
  implicit none
  integer(8), intent(in) :: backend, np, leafsize
  real(8),    intent(in) :: pts(3, np)      ! column = point (interleaved xyz)
  call kdtree_build(backend, pts, np, leafsize)
end subroutine axk_kdtree_build_r64

subroutine axk_kdtree_ball_r64(qpt, radius, nmax, out, nout)
  use kdtree_mod, only: kdtree_ball
  implicit none
  integer(8), intent(in)    :: nmax
  real(8),    intent(in)    :: qpt(3), radius
  integer(8), intent(inout) :: out(nmax)    ! 0-based original indices (filled up to nout)
  integer(8), intent(inout) :: nout         ! true count of points within radius
  call kdtree_ball(qpt, radius, out, nout)
end subroutine axk_kdtree_ball_r64

subroutine axk_kdtree_free_r64()
  use kdtree_mod, only: kdtree_free
  implicit none
  call kdtree_free()
end subroutine axk_kdtree_free_r64
