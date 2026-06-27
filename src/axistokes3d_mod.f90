module axistokes3d_mod
  ! ------------------------------------------------------------------
  ! Primitives module for AxiStokes3D: owns the package kind parameters
  ! (r64, r128, c64, c128) plus small generic utilities (gauss).  Other
  ! modules do
  !   use axistokes3d_mod, only: r64, r128, ...
  ! rather than re-defining the kinds themselves.
  ! ------------------------------------------------------------------
  implicit none

  integer, parameter :: r64  = 8
  integer, parameter :: r128 = 16
  integer, parameter :: c64  = 8
  integer, parameter :: c128 = 16

contains

  ! ------------------------------------------------------------------
  ! gauss
  ! Gauss-Legendre nodes tgl, weights wgl, and differentiation matrix Dgl
  ! on [-1,1].  Newton-Raphson on P_n (cosine initial guess), D-matrix via
  ! Fornberg formula.  Mirrors gauss.m (von Winckel / Bowei Wu).
  ! ------------------------------------------------------------------
  subroutine gauss_r64(n, tgl, wgl, Dgl)
    integer(8), intent(in)    :: n
    real(r64),  intent(inout) :: tgl(n), wgl(n), Dgl(n,n)

    integer(8) :: i, j, k, m
    real(r64)  :: z, z1, p0, p1, p2, pp, pi, tol
    real(r64)  :: a(n)

    pi  = acos(-1.0_r64)
    tol = epsilon(1.0_r64)
    m   = (n + 1_8) / 2_8

    ! ---- nodes and weights (Newton-Raphson, exploit symmetry) ----
    do i = 1, m
      z = cos(pi * (real(i,r64) - 0.25_r64) / (real(n,r64) + 0.5_r64))
      do
        p0 = 1.0_r64;  p1 = z
        do k = 2, n
          p2 = ((2.0_r64*real(k,r64) - 1.0_r64)*z*p1 &
               - (real(k,r64) - 1.0_r64)*p0) / real(k,r64)
          p0 = p1;  p1 = p2
        end do
        pp = real(n,r64) * (z*p1 - p0) / (z*z - 1.0_r64)
        z1 = z
        z  = z1 - p1/pp
        if (abs(z - z1) <= tol * max(1.0_r64, abs(z))) exit
      end do
      tgl(i)       = -z
      tgl(n+1-i)   =  z
      wgl(i)       = 2.0_r64 / ((1.0_r64 - z*z) * pp*pp)
      wgl(n+1-i)   = wgl(i)
    end do

    ! ---- differentiation matrix (Fornberg) ----
    do k = 1, n
      a(k) = 1.0_r64
      do j = 1, n
        if (j /= k) a(k) = a(k) * (tgl(k) - tgl(j))
      end do
    end do
    do k = 1, n
      do j = 1, n
        if (j /= k) then
          Dgl(j,k) = (a(j)/a(k)) / (tgl(j) - tgl(k))
        end if
      end do
      Dgl(k,k) = 0.0_r64
      do j = 1, n
        if (j /= k) Dgl(k,k) = Dgl(k,k) + 1.0_r64/(tgl(k) - tgl(j))
      end do
    end do

  end subroutine gauss_r64

  subroutine gauss_r128(n, tgl, wgl, Dgl)
    integer(8), intent(in)    :: n
    real(r128), intent(inout) :: tgl(n), wgl(n), Dgl(n,n)

    integer(8) :: i, j, k, m
    real(r128) :: z, z1, p0, p1, p2, pp, pi, tol
    real(r128) :: a(n)

    pi  = acos(-1.0_r128)
    tol = epsilon(1.0_r128)
    m   = (n + 1_8) / 2_8

    ! ---- nodes and weights ----
    do i = 1, m
      z = cos(pi * (real(i,r128) - 0.25_r128) / (real(n,r128) + 0.5_r128))
      do
        p0 = 1.0_r128;  p1 = z
        do k = 2, n
          p2 = ((2.0_r128*real(k,r128) - 1.0_r128)*z*p1 &
               - (real(k,r128) - 1.0_r128)*p0) / real(k,r128)
          p0 = p1;  p1 = p2
        end do
        pp = real(n,r128) * (z*p1 - p0) / (z*z - 1.0_r128)
        z1 = z
        z  = z1 - p1/pp
        if (abs(z - z1) <= tol * max(1.0_r128, abs(z))) exit
      end do
      tgl(i)       = -z
      tgl(n+1-i)   =  z
      wgl(i)       = 2.0_r128 / ((1.0_r128 - z*z) * pp*pp)
      wgl(n+1-i)   = wgl(i)
    end do

    ! ---- differentiation matrix (Fornberg) ----
    do k = 1, n
      a(k) = 1.0_r128
      do j = 1, n
        if (j /= k) a(k) = a(k) * (tgl(k) - tgl(j))
      end do
    end do
    do k = 1, n
      do j = 1, n
        if (j /= k) then
          Dgl(j,k) = (a(j)/a(k)) / (tgl(j) - tgl(k))
        end if
      end do
      Dgl(k,k) = 0.0_r128
      do j = 1, n
        if (j /= k) Dgl(k,k) = Dgl(k,k) + 1.0_r128/(tgl(k) - tgl(j))
      end do
    end do

  end subroutine gauss_r128

  ! ------------------------------------------------------------------
  ! ellipke_mc
  ! Complete elliptic K(m), E(m) and K-E via AGM, in the complementary
  ! parameter mc = 1 - m (cancellation-free as mc -> 0).  Shared by the
  ! axisym close-eval coef builders (directly or through carrier).
  ! ------------------------------------------------------------------
  subroutine ellipke_mc_r64(mc, K, E, KmE)
    real(r64), intent(in)    :: mc
    real(r64), intent(inout) :: K, E, KmE
    real(r64)  :: a, b, c0, s, p2, an, bn, cn, pi
    integer(8) :: it
    pi = acos(-1.0_r64)
    a = 1.0_r64; b = sqrt(mc); c0 = sqrt(1.0_r64 - mc); s = 0.5_r64*c0*c0; p2 = 1.0_r64
    do it = 1, 60
      an = (a + b)/2.0_r64; bn = sqrt(a*b); cn = (a - b)/2.0_r64
      a = an; b = bn; s = s + p2*cn*cn; p2 = p2*2.0_r64
      if (abs(cn) < 1.0e-17_r64) exit
    end do
    K = pi/(2.0_r64*a); E = K*(1.0_r64 - s); KmE = K*s
  end subroutine ellipke_mc_r64


  ! ------------------------------------------------------------------
  ! lagrange_interp
  ! 1D Lagrange interpolation matrix M(nt,ns): values at the source nodes
  ! xs(ns) -> values at the target nodes xt(nt), M(i,j) = prod_{k/=j} (xt_i-xs_k)/(xs_j-xs_k).
  ! Shared by the block builders for panel refinement, the p->2p close-eval
  ! upsample, and the refined->original fold-back projection.
  ! ------------------------------------------------------------------
  subroutine lagrange_interp_r64(ns, xs, nt, xt, M)
    integer(8), intent(in)    :: ns, nt
    real(r64),  intent(in)    :: xs(ns), xt(nt)
    real(r64),  intent(inout) :: M(nt,ns)
    integer(8) :: i, jj, kk
    real(r64)  :: num, den
    do i = 1, nt
      do jj = 1, ns
        num = 1.0_r64; den = 1.0_r64
        do kk = 1, ns
          if (kk /= jj) then
            num = num * (xt(i) - xs(kk))
            den = den * (xs(jj) - xs(kk))
          end if
        end do
        M(i,jj) = num/den
      end do
    end do
  end subroutine lagrange_interp_r64

end module axistokes3d_mod
