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

  ! ------------------------------------------------------------------
  ! sto3dslp_eval  (matrix-free 3D Stokes single-layer velocity)
  ! Same kernel as Sto3dSLPmat (utils/Sto3dSLPmat.m), without forming the
  ! dense 3M x 3N matrix:
  !   u_i(x) = 1/(8 pi) sum_s w_s [ f_i/r + (f.d) d_i / r^3 ],  d = x - y_s.
  ! tx(3,nt) targets, sx(3,ns) sources, sw(ns) weights, f(ns,3)=[fx fy fz] density,
  ! u(3,nt) velocity.  (mu = 1, matching Sto3dSLPmat.)
  ! ------------------------------------------------------------------
  subroutine sto3dslp_eval_r64(nt, tx, ns, sx, sw, f, u)
    integer(8), intent(in)    :: nt, ns
    real(r64),  intent(in)    :: tx(3,nt), sx(3,ns), sw(ns), f(ns,3)
    real(r64),  intent(inout) :: u(3,nt)
    integer(8) :: it, is
    real(r64)  :: d1, d2, d3, r2, ir, ir3, fd, w, ux, uy, uz, c

    c = 1.0_r64 / (8.0_r64 * acos(-1.0_r64))
    !$omp parallel do default(shared) schedule(static) &
    !$omp   private(it, is, d1, d2, d3, r2, ir, ir3, fd, w, ux, uy, uz)
    do it = 1, nt
      ux = 0.0_r64; uy = 0.0_r64; uz = 0.0_r64
      do is = 1, ns
        d1 = tx(1,it) - sx(1,is); d2 = tx(2,it) - sx(2,is); d3 = tx(3,it) - sx(3,is)
        r2 = d1*d1 + d2*d2 + d3*d3
        ir = 1.0_r64 / sqrt(r2);  ir3 = ir / r2
        w  = sw(is)
        fd = d1*f(is,1) + d2*f(is,2) + d3*f(is,3)
        ux = ux + w*(f(is,1)*ir + fd*d1*ir3)
        uy = uy + w*(f(is,2)*ir + fd*d2*ir3)
        uz = uz + w*(f(is,3)*ir + fd*d3*ir3)
      end do
      u(1,it) = c*ux;  u(2,it) = c*uy;  u(3,it) = c*uz
    end do
    !$omp end parallel do
  end subroutine sto3dslp_eval_r64

  ! ------------------------------------------------------------------
  ! sto3ddlp_eval  (matrix-free 3D Stokes double-layer / stresslet velocity)
  ! Same kernel as Sto3dDLPmat (utils/Sto3dDLPmat.m), no dense matrix:
  !   u_i(x) = 3/(4 pi) sum_s w_s (d.n)(d.f) d_i / r^5,  d = x - y_s, n = source normal.
  ! tx(3,nt) targets, sx(3,ns) sources, snx(3,ns) source normals, sw(ns) weights,
  ! f(ns,3)=[fx fy fz] density, u(3,nt) velocity.  (mu-independent.)
  ! ------------------------------------------------------------------
  subroutine sto3ddlp_eval_r64(nt, tx, ns, sx, snx, sw, f, u)
    integer(8), intent(in)    :: nt, ns
    real(r64),  intent(in)    :: tx(3,nt), sx(3,ns), snx(3,ns), sw(ns), f(ns,3)
    real(r64),  intent(inout) :: u(3,nt)
    integer(8) :: it, is
    real(r64)  :: d1, d2, d3, r2, ir, ir5, dn, df, g, ux, uy, uz, c

    c = 3.0_r64 / (4.0_r64 * acos(-1.0_r64))
    !$omp parallel do default(shared) schedule(static) &
    !$omp   private(it, is, d1, d2, d3, r2, ir, ir5, dn, df, g, ux, uy, uz)
    do it = 1, nt
      ux = 0.0_r64; uy = 0.0_r64; uz = 0.0_r64
      do is = 1, ns
        d1 = tx(1,it) - sx(1,is); d2 = tx(2,it) - sx(2,is); d3 = tx(3,it) - sx(3,is)
        r2 = d1*d1 + d2*d2 + d3*d3
        ir = 1.0_r64 / sqrt(r2);  ir5 = ir / (r2*r2)
        dn = d1*snx(1,is) + d2*snx(2,is) + d3*snx(3,is)
        df = d1*f(is,1) + d2*f(is,2) + d3*f(is,3)
        g  = sw(is) * dn * df * ir5
        ux = ux + g*d1;  uy = uy + g*d2;  uz = uz + g*d3
      end do
      u(1,it) = c*ux;  u(2,it) = c*uy;  u(3,it) = c*uz
    end do
    !$omp end parallel do
  end subroutine sto3ddlp_eval_r64

end module axistokes3d_mod
