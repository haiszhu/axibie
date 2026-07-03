! kdtree_mod.f90
!
! Fortran front-end for C++ kd-tree libraries (bridge in kdtree_cwrap.cpp).
! Stateful build-once / query-many sequence (one tree held at a time):
!     call kdtree_build(backend, pts, npts, leafsize)   ! backend = KD_NBODY | KD_TAIYA (required)
!     call kdtree_ball (qpt, radius, out, nout)          ! all points within radius of qpt
!     call kdtree_free ()
!   nbodyhpc : built in C++, exported to module arrays, radius traversal done HERE in Fortran.
!   taiya    : built + queried in C++ (native ball query); tree retained in the C++ bridge.
! `out` returns 0-based ORIGINAL point indices, nout = count (true count; may exceed size(out)).
!
! Boundary convention (-i8): integer*8 <-> int64_t, real*4 <-> float, real*8 <-> double.

module kdtree_mod
  implicit none
  private

  public :: kdtree_build, kdtree_ball, kdtree_free
  public :: KD_NBODY, KD_TAIYA

  integer*8, parameter :: KD_NBODY = 1
  integer*8, parameter :: KD_TAIYA = 2

  ! ---- module state: the currently-held tree ----
  integer*8, save :: s_backend = 0
  integer*8, allocatable, save :: s_dim(:), s_lft(:), s_rgt(:), s_idx(:)  ! nbody exported tree
  real*4,    allocatable, save :: s_split(:), s_xp(:), s_yp(:), s_zp(:)
  !$omp threadprivate(s_backend, s_dim, s_lft, s_rgt, s_idx, s_split, s_xp, s_yp, s_zp)

  ! ---- C bridge (kdtree_cwrap.cpp; external symbols, gfortran appends "_") ----
  interface
    subroutine kdtree_nbody_build(np, pts, leafsize, maxthreads, nnodes, npad)
      implicit none
      integer*8, intent(in)  :: np, leafsize, maxthreads
      real*4,    intent(in)  :: pts(*)
      integer*8, intent(out) :: nnodes, npad
    end subroutine
    subroutine kdtree_nbody_export(dim, split, lft, rgt, xp, yp, zp, idx)
      implicit none
      integer*8, intent(out) :: dim(*), lft(*), rgt(*), idx(*)
      real*4,    intent(out) :: split(*), xp(*), yp(*), zp(*)
    end subroutine
    subroutine kdtree_nbody_free()
    end subroutine

    subroutine kdtree_taiya_build(np, pts)
      implicit none
      integer*8, intent(in) :: np
      real*8,    intent(in) :: pts(*)
    end subroutine
    subroutine kdtree_taiya_ball(qpt, radius, nmax, out, nout)
      implicit none
      integer*8, intent(in)  :: nmax
      real*8,    intent(in)  :: qpt(3), radius
      integer*8, intent(out) :: out(*), nout
    end subroutine
    subroutine kdtree_taiya_free()
    end subroutine
  end interface

contains

  ! Build the tree from pts(3,npts) (column = point) with the chosen backend (required).
  subroutine kdtree_build(backend, pts, npts, leafsize)
    implicit none
    integer*8, intent(in) :: backend, npts, leafsize
    real*8,    intent(in) :: pts(3, npts)
    integer*8 :: nnodes, npad, mt, i
    real*4, allocatable :: ptsf(:)
    real*8, allocatable :: ptsd(:)

    call kdtree_free()          ! drop any previously-held tree
    s_backend = backend
    mt = 1_8

    if (backend == KD_NBODY) then
      allocate(ptsf(3*npts))
      do i = 1, npts
        ptsf(3*i-2) = real(pts(1,i), 4); ptsf(3*i-1) = real(pts(2,i), 4); ptsf(3*i) = real(pts(3,i), 4)
      end do
      call kdtree_nbody_build(npts, ptsf, leafsize, mt, nnodes, npad)
      allocate(s_dim(nnodes), s_split(nnodes), s_lft(nnodes), s_rgt(nnodes))
      allocate(s_xp(npad), s_yp(npad), s_zp(npad), s_idx(npad))
      call kdtree_nbody_export(s_dim, s_split, s_lft, s_rgt, s_xp, s_yp, s_zp, s_idx)
      call kdtree_nbody_free()  ! exported copy is all we need
      deallocate(ptsf)

    else if (backend == KD_TAIYA) then
      allocate(ptsd(3*npts))
      do i = 1, npts
        ptsd(3*i-2) = pts(1,i); ptsd(3*i-1) = pts(2,i); ptsd(3*i) = pts(3,i)
      end do
      call kdtree_taiya_build(npts, ptsd)
      deallocate(ptsd)
    end if
  end subroutine kdtree_build

  ! All points within `radius` of `qpt(3)`: 0-based original indices in out(1:nout).
  subroutine kdtree_ball(qpt, radius, out, nout)
    implicit none
    real*8,    intent(in)  :: qpt(3), radius
    integer*8, intent(out) :: out(:), nout
    integer*8 :: nmax
    real*4    :: q4(3), r4

    nmax = size(out, kind=8)
    nout = 0_8

    if (s_backend == KD_NBODY) then
      q4 = real(qpt, 4); r4 = real(radius, 4)
      call ball_bbox_query(0_8, s_dim, s_split, s_lft, s_rgt, s_xp, s_yp, s_zp, s_idx, &
                           q4, r4, r4*r4, out, nout, nmax)
    else if (s_backend == KD_TAIYA) then
      call kdtree_taiya_ball(qpt, radius, nmax, out, nout)
    end if
  end subroutine kdtree_ball

  ! Release the held tree.
  subroutine kdtree_free()
    implicit none
    if (s_backend == KD_TAIYA) call kdtree_taiya_free()
    if (allocated(s_dim))   deallocate(s_dim, s_split, s_lft, s_rgt)
    if (allocated(s_xp))    deallocate(s_xp, s_yp, s_zp, s_idx)
    s_backend = 0_8
  end subroutine kdtree_free

  ! Recursive radius traversal over the exported flat nbodyhpc tree (0-based node index).
  recursive subroutine ball_bbox_query(node, dim, split, lft, rgt, &
                                       xp, yp, zp, idx, q, r, r2, out, nout, nmax)
    implicit none
    integer*8, intent(in)    :: node, nmax
    integer*8, intent(in)    :: dim(:), lft(:), rgt(:), idx(:)
    real*4,    intent(in)    :: split(:), xp(:), yp(:), zp(:)
    real*4,    intent(in)    :: q(3), r, r2
    integer*8, intent(inout) :: out(:), nout
    integer*8 :: s, nd
    real*4    :: dx, dy, dz, sp

    nd = dim(node + 1)
    if (nd == -1_8) then
      do s = lft(node + 1), rgt(node + 1) - 1_8
        dx = xp(s+1) - q(1); dy = yp(s+1) - q(2); dz = zp(s+1) - q(3)
        if (dx*dx + dy*dy + dz*dz <= r2) then
          nout = nout + 1_8
          if (nout <= nmax) out(nout) = idx(s+1)
        end if
      end do
    else
      sp = split(node + 1)
      if (q(nd+1) - r <= sp) &
        call ball_bbox_query(lft(node+1), dim,split,lft,rgt, xp,yp,zp,idx, q,r,r2, out,nout,nmax)
      if (q(nd+1) + r >= sp) &
        call ball_bbox_query(rgt(node+1), dim,split,lft,rgt, xp,yp,zp,idx, q,r,r2, out,nout,nmax)
    end if
  end subroutine ball_bbox_query

end module kdtree_mod
