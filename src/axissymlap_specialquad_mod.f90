module axissymlap_specialquad_mod
  ! ==================================================================
  ! Per-panel 0th-mode close-eval OPERATOR for the axisymmetric Laplace ladder
  ! S, S', D, D', D'' (SLP/SLPn/DLP/DLPn/DLPnn).  Laplace analog of
  ! axissymstok_specialquad_mod / the axm_specialquad_* family: assembles the
  ! nt x p close-eval block W by combining the 0th-mode kernel-split coefficients
  ! (axissymlap_kernelsplit_mod) with the Helsing singular quadrature
  ! (specialquad_mod: As log, Ad Cauchy, Az=A1-iA2 deriv-Cauchy, Azz=A3-iA4 cubic)
  ! and the caller's panel quadrature (q->p Legendre projection Pmat, upsampled
  ! weights gw).  Scalar (nt x p, not 3nt x 3p) and 0th-mode (no mode loop).
  !
  ! Geometry pre-evaluated by the host (Fortran cannot take MATLAB handles):
  ! upsampled source nodes y = Z(tt) (q=2p), tangents yp = Zp(tt), endpoints
  ! za=Z(a0), zb=Z(b0), parameter half-length hm=(b0-a0)/2.  The source normal
  ! nv = -i yp/|yp| is derived here.  Coefficients are REAL; the normals nu=nx
  ! (target) and nu'=nv (source) are folded into the special-quad weights exactly
  ! as the MATLAB Lap*Axi0Specialquad_v2 references do:
  !   target Cauchy   Re[nu  Ad conj(nv)]     (bare Cauchy = Ad/nv = Ad conj(nv))
  !   source Cauchy   Re[Ad]                   (nu' built into the DLP value Ad)
  !   nu nu' deriv-C  Re[nu  Az]               (nu' built into Az)
  !   nu^2  deriv-C   Re[-nu^2 conj(nv) Az]    (conj(nv) divides nu' out)
  !   cubic           2 Re[nu^2 Azz]           (nu' built into Azz)
  !
  !   W(nt,p) = blk(nt,q) * Pmat(q,p)
  ! ==================================================================
  use axilaplace3d_mod, only: r64
  use axissymlap_kernelsplit_mod, only: &
       coef_slp0   => axissymlap_kernelsplit_coef_slp0_r64,   &
       coef_slpn0  => axissymlap_kernelsplit_coef_slpn0_r64,  &
       coef_dlp0   => axissymlap_kernelsplit_coef_dlp0_r64,   &
       coef_dlpn0  => axissymlap_kernelsplit_coef_dlpn0_r64,  &
       coef_dlpnn0 => axissymlap_kernelsplit_coef_dlpnn0_r64
  use specialquad_mod, only: sdspecialquad_r64
  implicit none
  private
  public :: axissymlap_specialquad_slp0_r64
  public :: axissymlap_specialquad_slpn0_r64
  public :: axissymlap_specialquad_dlp0_r64
  public :: axissymlap_specialquad_dlpn0_r64
  public :: axissymlap_specialquad_dlpnn0_r64

contains

  ! ---- SLP : log + smooth ----
  subroutine axissymlap_specialquad_slp0_r64(nt, z, p, hm, y, yp, za, zb, gw, Pmat, iside, W)
    integer(8),   intent(in)    :: nt, p, iside
    real(r64),    intent(in)    :: hm, gw(2*p), Pmat(2*p,p)
    complex(r64), intent(in)    :: z(nt), y(2*p), yp(2*p), za, zb
    real(r64),    intent(inout) :: W(nt,p)
    integer(8) :: q, i, j
    real(r64)  :: pi, twopi, spdj
    complex(r64) :: ci
    real(r64),    allocatable :: C1(:,:), C2(:,:), As(:,:), A1(:,:), A2(:,:), A3(:,:), A4(:,:), arc(:), blk(:,:)
    complex(r64), allocatable :: Ad(:,:), nv(:), wxp(:)
    pi = acos(-1.0_r64); twopi = 2.0_r64*pi; ci = (0.0_r64,1.0_r64); q = 2*p
    allocate(C1(nt,q),C2(nt,q),As(q,nt),Ad(q,nt),A1(q,nt),A2(q,nt),A3(q,nt),A4(q,nt))
    allocate(arc(q),nv(q),wxp(q),blk(nt,q))
    do j = 1, q
      spdj = abs(yp(j)); nv(j) = -ci*yp(j)/spdj; arc(j) = spdj*gw(j)*hm; wxp(j) = yp(j)*gw(j)*hm
    end do
    call coef_slp0(nt, q, real(z,r64),aimag(z), real(y,r64),aimag(y), C1, C2)
    call sdspecialquad_r64(nt, z, q, y, nv, wxp, za, zb, iside, As, Ad, A1, A2, A3, A4)
    do i = 1, nt
      do j = 1, q
        blk(i,j) = twopi*C1(i,j)*As(j,i) + C2(i,j)*arc(j)
      end do
    end do
    W = matmul(blk, Pmat)
    deallocate(C1,C2,As,Ad,A1,A2,A3,A4,arc,nv,wxp,blk)
  end subroutine axissymlap_specialquad_slp0_r64

  ! ---- SLPn : + target Cauchy.  nx = target normal (nt). ----
  subroutine axissymlap_specialquad_slpn0_r64(nt, z, nx, p, hm, y, yp, za, zb, gw, Pmat, iside, W)
    integer(8),   intent(in)    :: nt, p, iside
    real(r64),    intent(in)    :: hm, gw(2*p), Pmat(2*p,p)
    complex(r64), intent(in)    :: z(nt), nx(nt), y(2*p), yp(2*p), za, zb
    real(r64),    intent(inout) :: W(nt,p)
    integer(8) :: q, i, j
    real(r64)  :: pi, twopi, spdj
    complex(r64) :: ci
    real(r64),    allocatable :: C1(:,:), C2(:,:), C3(:,:), As(:,:), A1(:,:), A2(:,:), A3(:,:), A4(:,:), arc(:), blk(:,:)
    complex(r64), allocatable :: Ad(:,:), nv(:), wxp(:)
    pi = acos(-1.0_r64); twopi = 2.0_r64*pi; ci = (0.0_r64,1.0_r64); q = 2*p
    allocate(C1(nt,q),C2(nt,q),C3(nt,q),As(q,nt),Ad(q,nt),A1(q,nt),A2(q,nt),A3(q,nt),A4(q,nt))
    allocate(arc(q),nv(q),wxp(q),blk(nt,q))
    do j = 1, q
      spdj = abs(yp(j)); nv(j) = -ci*yp(j)/spdj; arc(j) = spdj*gw(j)*hm; wxp(j) = yp(j)*gw(j)*hm
    end do
    call coef_slpn0(nt, q, real(z,r64),aimag(z), real(y,r64),aimag(y), real(nx,r64),aimag(nx), C1, C2, C3)
    call sdspecialquad_r64(nt, z, q, y, nv, wxp, za, zb, iside, As, Ad, A1, A2, A3, A4)
    do i = 1, nt
      do j = 1, q
        blk(i,j) = twopi*C1(i,j)*As(j,i) + C2(i,j)*arc(j) &
                 - twopi*C3(i,j)*real(nx(i)*Ad(j,i)*conjg(nv(j)), r64)   ! Az_S = -Ad conj(nv)
      end do
    end do
    W = matmul(blk, Pmat)
    deallocate(C1,C2,C3,As,Ad,A1,A2,A3,A4,arc,nv,wxp,blk)
  end subroutine axissymlap_specialquad_slpn0_r64

  ! ---- DLP : + source Cauchy (nu' built into Ad via the geometric source normal nv). ----
  subroutine axissymlap_specialquad_dlp0_r64(nt, z, p, hm, y, yp, za, zb, gw, Pmat, iside, W)
    integer(8),   intent(in)    :: nt, p, iside
    real(r64),    intent(in)    :: hm, gw(2*p), Pmat(2*p,p)
    complex(r64), intent(in)    :: z(nt), y(2*p), yp(2*p), za, zb
    real(r64),    intent(inout) :: W(nt,p)
    integer(8) :: q, i, j
    real(r64)  :: pi, twopi, spdj
    complex(r64) :: ci
    real(r64),    allocatable :: C1(:,:), C2(:,:), C3(:,:), As(:,:), A1(:,:), A2(:,:), A3(:,:), A4(:,:), arc(:), blk(:,:)
    complex(r64), allocatable :: Ad(:,:), nv(:), wxp(:)
    pi = acos(-1.0_r64); twopi = 2.0_r64*pi; ci = (0.0_r64,1.0_r64); q = 2*p
    allocate(C1(nt,q),C2(nt,q),C3(nt,q),As(q,nt),Ad(q,nt),A1(q,nt),A2(q,nt),A3(q,nt),A4(q,nt))
    allocate(arc(q),nv(q),wxp(q),blk(nt,q))
    do j = 1, q
      spdj = abs(yp(j)); nv(j) = -ci*yp(j)/spdj; arc(j) = spdj*gw(j)*hm; wxp(j) = yp(j)*gw(j)*hm
    end do
    call coef_dlp0(nt, q, real(z,r64),aimag(z), real(y,r64),aimag(y), real(nv,r64),aimag(nv), C1, C2, C3)
    call sdspecialquad_r64(nt, z, q, y, nv, wxp, za, zb, iside, As, Ad, A1, A2, A3, A4)
    do i = 1, nt
      do j = 1, q
        blk(i,j) = twopi*C1(i,j)*As(j,i) + C2(i,j)*arc(j) &
                 + twopi*C3(i,j)*real(Ad(j,i), r64)
      end do
    end do
    W = matmul(blk, Pmat)
    deallocate(C1,C2,C3,As,Ad,A1,A2,A3,A4,arc,nv,wxp,blk)
  end subroutine axissymlap_specialquad_dlp0_r64

  ! ---- DLPn : target Cauchy + source Cauchy + (nu nu') deriv-Cauchy. ----
  subroutine axissymlap_specialquad_dlpn0_r64(nt, z, nx, p, hm, y, yp, za, zb, gw, Pmat, iside, W)
    integer(8),   intent(in)    :: nt, p, iside
    real(r64),    intent(in)    :: hm, gw(2*p), Pmat(2*p,p)
    complex(r64), intent(in)    :: z(nt), nx(nt), y(2*p), yp(2*p), za, zb
    real(r64),    intent(inout) :: W(nt,p)
    integer(8) :: q, i, j
    real(r64)  :: pi, twopi, spdj
    complex(r64) :: ci, Az
    real(r64),    allocatable :: C1(:,:), C2(:,:), C3a(:,:), C3b(:,:), C4(:,:)
    real(r64),    allocatable :: As(:,:), A1(:,:), A2(:,:), A3(:,:), A4(:,:), arc(:), blk(:,:)
    complex(r64), allocatable :: Ad(:,:), nv(:), wxp(:)
    pi = acos(-1.0_r64); twopi = 2.0_r64*pi; ci = (0.0_r64,1.0_r64); q = 2*p
    allocate(C1(nt,q),C2(nt,q),C3a(nt,q),C3b(nt,q),C4(nt,q))
    allocate(As(q,nt),Ad(q,nt),A1(q,nt),A2(q,nt),A3(q,nt),A4(q,nt),arc(q),nv(q),wxp(q),blk(nt,q))
    do j = 1, q
      spdj = abs(yp(j)); nv(j) = -ci*yp(j)/spdj; arc(j) = spdj*gw(j)*hm; wxp(j) = yp(j)*gw(j)*hm
    end do
    call coef_dlpn0(nt, q, real(z,r64),aimag(z), real(y,r64),aimag(y), &
         real(nx,r64),aimag(nx), real(nv,r64),aimag(nv), C1, C2, C3a, C3b, C4)
    call sdspecialquad_r64(nt, z, q, y, nv, wxp, za, zb, iside, As, Ad, A1, A2, A3, A4)
    do i = 1, nt
      do j = 1, q
        Az = cmplx(A1(j,i), -A2(j,i), r64)
        blk(i,j) = twopi*C1(i,j)*As(j,i) + C2(i,j)*arc(j) &
                 - twopi*C3a(i,j)*real(nx(i)*Ad(j,i)*conjg(nv(j)), r64) &   ! Az_S = -Ad conj(nv)
                 + twopi*C3b(i,j)*real(Ad(j,i), r64) &
                 + twopi*C4(i,j)*real(nx(i)*Az, r64)
      end do
    end do
    W = matmul(blk, Pmat)
    deallocate(C1,C2,C3a,C3b,C4,As,Ad,A1,A2,A3,A4,arc,nv,wxp,blk)
  end subroutine axissymlap_specialquad_dlpn0_r64

  ! ---- DLPnn : target/source Cauchy + (nu^2) & (nu nu') deriv-Cauchy + cubic. ----
  subroutine axissymlap_specialquad_dlpnn0_r64(nt, z, nx, p, hm, y, yp, za, zb, gw, Pmat, iside, W)
    integer(8),   intent(in)    :: nt, p, iside
    real(r64),    intent(in)    :: hm, gw(2*p), Pmat(2*p,p)
    complex(r64), intent(in)    :: z(nt), nx(nt), y(2*p), yp(2*p), za, zb
    real(r64),    intent(inout) :: W(nt,p)
    integer(8) :: q, i, j
    real(r64)  :: pi, twopi, fourpi, spdj
    complex(r64) :: ci, Az, Azz, nxi, nxi2
    real(r64),    allocatable :: C1(:,:), C2(:,:), C3a(:,:), C3b(:,:), C4a(:,:), C4b(:,:), C5(:,:)
    real(r64),    allocatable :: As(:,:), A1(:,:), A2(:,:), A3(:,:), A4(:,:), arc(:), blk(:,:)
    complex(r64), allocatable :: Ad(:,:), nv(:), wxp(:)
    pi = acos(-1.0_r64); twopi = 2.0_r64*pi; fourpi = 4.0_r64*pi; ci = (0.0_r64,1.0_r64); q = 2*p
    allocate(C1(nt,q),C2(nt,q),C3a(nt,q),C3b(nt,q),C4a(nt,q),C4b(nt,q),C5(nt,q))
    allocate(As(q,nt),Ad(q,nt),A1(q,nt),A2(q,nt),A3(q,nt),A4(q,nt),arc(q),nv(q),wxp(q),blk(nt,q))
    do j = 1, q
      spdj = abs(yp(j)); nv(j) = -ci*yp(j)/spdj; arc(j) = spdj*gw(j)*hm; wxp(j) = yp(j)*gw(j)*hm
    end do
    call coef_dlpnn0(nt, q, real(z,r64),aimag(z), real(y,r64),aimag(y), &
         real(nx,r64),aimag(nx), real(nv,r64),aimag(nv), C1, C2, C3a, C3b, C4a, C4b, C5)
    call sdspecialquad_r64(nt, z, q, y, nv, wxp, za, zb, iside, As, Ad, A1, A2, A3, A4)
    do i = 1, nt
      nxi = nx(i); nxi2 = nxi*nxi
      do j = 1, q
        Az  = cmplx(A1(j,i), -A2(j,i), r64)
        Azz = cmplx(A3(j,i), -A4(j,i), r64)
        blk(i,j) = twopi*C1(i,j)*As(j,i) + C2(i,j)*arc(j) &
                 - fourpi*C3a(i,j)*real(nxi*Ad(j,i)*conjg(nv(j)), r64) &    ! Az_S = -Ad conj(nv)
                 + twopi *C3b(i,j)*real(Ad(j,i), r64) &
                 + twopi *C4a(i,j)*real(-nxi2*conjg(nv(j))*Az, r64) &
                 + fourpi*C4b(i,j)*real(nxi*Az, r64) &
                 + fourpi*C5(i,j) *real(nxi2*Azz, r64)
      end do
    end do
    W = matmul(blk, Pmat)
    deallocate(C1,C2,C3a,C3b,C4a,C4b,C5,As,Ad,A1,A2,A3,A4,arc,nv,wxp,blk)
  end subroutine axissymlap_specialquad_dlpnn0_r64

end module axissymlap_specialquad_mod
