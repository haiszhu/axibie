! specialquad_mod.f90
! Helsing close-evaluation singular quadrature for 2D Laplace BIE on panels.
! Port of SDspecialquad from GRF_disk.m / BIE2D (Barnett, Helsing 2009).
! Singular-quadrature companion to the kernel-split coefficients: As=Sspecialquad value (log/SLP,
! pairs with C1), Ad=Dspecialquad value (Cauchy, C3), Az=A1/A2 (deriv-Cauchy, C4),
! Azz=A3/A4 (deriv^2-Cauchy, C5/traction).
!
! For nt near targets and one panel [za,zb] with p source nodes, returns:
!   As(p,nt)  : SLP correction weights (real)
!   Ad(p,nt)   : DLP correction weights (complex; take real for Laplace)
!   A1(p,nt)  : Re(Az) = x-derivative DLP weights (real)
!   A2(p,nt)  : -Im(Az) = y-derivative DLP weights (real)
!   A3(p,nt)  : Re(Azz) (real)
!   A4(p,nt)  : -Im(Azz) (real)
!
! Note: Fortran outputs are p×nt (column-major).
!       MATLAB's SDspecialquad returns nt×p (row-major convention).
!       Transpose in MATLAB: As_matlab = As_fortran' (for nt=1, shape is the same).
!
! Usage in Laplace BIE for near panel kk, targets m=1..nt:
!   u(m) += dot(ub_kk, Ad(:,m)) - dot(unb_kk, As(:,m))
!   (use real(Ad(:,m)) for purely real Laplace problem)
!
! Internal helpers (private):
!   zge_r64 (p, Ad, zb) — p×p complex GE with partial pivoting, in-place
!   zge_r128(p, Ad, zb)
!
! Public routines:
!   sdspecialquad_r64 (nt, zt, p, zsrc, nzsrc, wzp, za, zb, iside, As, Ad, A1, A2, A3, A4)
!   sdspecialquad_r128(nt, zt, p, zsrc, nzsrc, wzp, za, zb, iside, As, Ad, A1, A2, A3, A4)

module specialquad_mod
  use axistokes3d_mod, only: r64, r128
  implicit none
  private
  public :: sdspecialquad_r64, sdspecialquad_r128

contains

  ! ------------------------------------------------------------------
  ! zge_r64: in-place LU with partial pivoting. Solves Ad*x = zb.
  ! On exit: zb = solution, Ad = destroyed.
  ! ------------------------------------------------------------------
  subroutine zge_r64(p, Ad, zb)
    integer(8),   intent(in)    :: p
    complex(r64), intent(inout) :: Ad(p,p), zb(p)
    integer(8)   :: i, j, k, ipiv
    complex(r64) :: tmp, fac

    do k = 1, p
      ipiv = k
      do i = k+1, p
        if (abs(Ad(i,k)) > abs(Ad(ipiv,k))) ipiv = i
      end do
      if (ipiv /= k) then
        tmp = zb(k);  zb(k) = zb(ipiv);  zb(ipiv) = tmp
        do j = 1, p
          tmp = Ad(k,j);  Ad(k,j) = Ad(ipiv,j);  Ad(ipiv,j) = tmp
        end do
      end if
      do i = k+1, p
        fac    = Ad(i,k) / Ad(k,k)
        Ad(i,k) = cmplx(0.0_r64, 0.0_r64, r64)
        do j = k+1, p
          Ad(i,j) = Ad(i,j) - fac*Ad(k,j)
        end do
        zb(i) = zb(i) - fac*zb(k)
      end do
    end do
    do k = p, 1, -1
      do j = k+1, p
        zb(k) = zb(k) - Ad(k,j)*zb(j)
      end do
      zb(k) = zb(k) / Ad(k,k)
    end do
  end subroutine zge_r64

  subroutine zge_r128(p, Ad, zb)
    integer(8),    intent(in)    :: p
    complex(r128), intent(inout) :: Ad(p,p), zb(p)
    integer(8)    :: i, j, k, ipiv
    complex(r128) :: tmp, fac

    do k = 1, p
      ipiv = k
      do i = k+1, p
        if (abs(Ad(i,k)) > abs(Ad(ipiv,k))) ipiv = i
      end do
      if (ipiv /= k) then
        tmp = zb(k);  zb(k) = zb(ipiv);  zb(ipiv) = tmp
        do j = 1, p
          tmp = Ad(k,j);  Ad(k,j) = Ad(ipiv,j);  Ad(ipiv,j) = tmp
        end do
      end if
      do i = k+1, p
        fac    = Ad(i,k) / Ad(k,k)
        Ad(i,k) = cmplx(0.0_r128, 0.0_r128, r128)
        do j = k+1, p
          Ad(i,j) = Ad(i,j) - fac*Ad(k,j)
        end do
        zb(i) = zb(i) - fac*zb(k)
      end do
    end do
    do k = p, 1, -1
      do j = k+1, p
        zb(k) = zb(k) - Ad(k,j)*zb(j)
      end do
      zb(k) = zb(k) / Ad(k,k)
    end do
  end subroutine zge_r128

  ! ------------------------------------------------------------------
  ! sdspecialquad_r64
  ! Close-evaluation correction weights for one panel [za,zb], nt targets.
  ! Port of SDspecialquad from GRF_disk.m (Barnett/Helsing 2009).
  !
  ! Arguments:
  !   nt         : number of targets
  !   zt(nt)     : target positions (complex r64)
  !   p          : panel order (# source nodes)
  !   zsrc(p)    : source nodes  (= s.x in MATLAB)
  !   nzsrc(p)   : source normals (= s.nx)
  !   wzp(p)     : complex speed weights w*zp (= s.wxp)
  !   za, zb       : panel endpoints (complex r64)
  !   iside      : 1 = exterior ('e'), 0 = interior ('i')
  !   As(p,nt)   : SLP weights (real r64)
  !   Ad(p,nt)    : DLP weights (complex r64; take real for Laplace)
  !   A1(p,nt)   : Re(Az), x-deriv of DLP (real r64)
  !   A2(p,nt)   : -Im(Az), y-deriv of DLP (real r64)
  !   A3(p,nt)   : Re(Azz) (real r64)
  !   A4(p,nt)   : -Im(Azz) (real r64)
  !
  ! Implementation note:
  !   Only the upward-recurrence path (|x| <= 100) is implemented.
  !   For list2 near-panel targets |x| << 1, this is always valid.
  !   The far-quadrature path (|x| > 100) is not needed for BIE list2.
  ! ------------------------------------------------------------------
  subroutine sdspecialquad_r64(nt, zt, p, zsrc, nzsrc, wzp, za, zb, iside, &
      As, Ad, A1, A2, A3, A4)

    integer(8),   intent(in)  :: nt, p
    complex(r64), intent(in)  :: zt(nt), zsrc(p), nzsrc(p), wzp(p), za, zb
    integer(8),   intent(in)  :: iside
    real(r64),    intent(out) :: As(p,nt), A1(p,nt), A2(p,nt), A3(p,nt), A4(p,nt)
    complex(r64), intent(out) :: Ad(p,nt)

    complex(r64) :: zsc, zmid, x, y(p), gam, log1mx, logprod
    complex(r64) :: c(p), V(p,p), VT(p,p)
    complex(r64) :: Pvec(p+1), Qvec(p), Rvec(p), Svec(p)
    complex(r64) :: cs(p), cd(p), cr(p), cs2(p)
    complex(r64) :: ci, Az_j, Azz_j
    real(r64)    :: pi, abszsc
    integer(8)   :: j, k, m

    pi = acos(-1.0_r64)
    ci = cmplx(0.0_r64, 1.0_r64, r64)

    zsc    = (zb - za) / 2.0_r64
    zmid   = (zb + za) / 2.0_r64
    abszsc = abs(zsc)

    ! transformed source nodes
    do j = 1, p
      y(j) = (zsrc(j) - zmid) / zsc
    end do

    ! Helsing c_k = (1 - (-1)^k)/k
    do k = 1, p
      if (mod(k,2) == 1) then
        c(k) = cmplx(2.0_r64 / real(k,r64), 0.0_r64, r64)
      else
        c(k) = cmplx(0.0_r64, 0.0_r64, r64)
      end if
    end do

    ! Vandermonde V(:,k) = y.^(k-1)
    do j = 1, p
      V(j,1) = cmplx(1.0_r64, 0.0_r64, r64)
      do k = 2, p
        V(j,k) = V(j,k-1) * y(j)
      end do
    end do

    ! MATLAB uses gam = exp(1i*pi/4); if side=='e', gam = conj(gam)
    gam = exp(ci * pi / 4.0_r64)
    if (iside == 1) gam = conjg(gam)

    do m = 1, nt

      x = (zt(m) - zmid) / zsc

      ! P(1) = log(gam) + log((1-x)/(gam*(-1-x)))
      log1mx  = log(gam) + log((1.0_r64 - x) / (gam * (-1.0_r64 - x)))
      Pvec(1) = log1mx

      if (abs(x) <= 1.75_r64) then
        ! NEAR: upward recurrence  P(k+1) = x*P(k) + c(k)
        do k = 1, p
          Pvec(k+1) = x * Pvec(k) + c(k)
        end do
      else
        ! FAR (|x|>1.75): seed top moment by direct quadrature, recurse DOWNWARD
        ! (Sspecialquad/Dspecialquad far branch).  P(p+1) = sum_j (w*y'/zsc) y_j^p/(y_j-x).
        Pvec(p+1) = cmplx(0.0_r64, 0.0_r64, r64)
        do j = 1, p
          Pvec(p+1) = Pvec(p+1) + (wzp(j)/zsc) * (V(j,p)*y(j)) / (y(j) - x)
        end do
        do k = p, 2, -1
          Pvec(k) = (Pvec(k+1) - c(k)) / x
        end do
      end if

      ! Q(odd)  = P(even) - log((1-x)*(-1-x))
      ! Q(even) = P(odd next) - [log(gam)+log((1-x)/(gam*(-1-x)))]
      logprod = log((1.0_r64 - x) * (-1.0_r64 - x))
      do k = 1, p
        if (mod(k,2) == 1) then
          Qvec(k) = (Pvec(k+1) - logprod) / real(k, r64)
        else
          Qvec(k) = (Pvec(k+1) - log1mx) / real(k, r64)
        end if
      end do

      ! solve V.' \ Qvec
      do j = 1, p
        do k = 1, p
          VT(k,j) = V(j,k)
        end do
      end do
      cs = Qvec
      call zge_r64(p, VT, cs)

      ! MATLAB:
      ! As = real((V.'\Q).'.*((1i*s.nx)')*zsc)/(2*pi*abs(zsc));
      ! As = As*abs(zsc) - log(abs(zsc))/(2*pi)*abs(s.wxp)';
      do j = 1, p
        ! Match MATLAB SDspecialquad exactly: (1i*s.nx)' is conjugate-transpose.
        As(j,m) = real(cs(j) * conjg(ci*nzsrc(j)) * zsc, r64) / (2.0_r64*pi*abszsc)
      end do
      do j = 1, p
        As(j,m) = As(j,m) * abszsc - log(abszsc)/(2.0_r64*pi) * abs(wzp(j))
      end do

      ! solve V.' \ P(1:p)
      do j = 1, p
        do k = 1, p
          VT(k,j) = V(j,k)
        end do
      end do
      cd = Pvec(1:p)
      call zge_r64(p, VT, cd)

      do j = 1, p
        Ad(j,m) = cd(j) * (ci / (2.0_r64*pi))
      end do

      ! R = -(1/(1-x) + (-1)^(k-1)/(1+x)) + (k-1)*[0;P(1:p-1)]
      do k = 1, p
        if (k == 1) then
          Rvec(k) = -( 1.0_r64/(1.0_r64 - x) + 1.0_r64/(1.0_r64 + x) )
        else
          if (mod(k-1,2) == 0) then
            Rvec(k) = -( 1.0_r64/(1.0_r64 - x) + 1.0_r64/(1.0_r64 + x) ) &
                      + real(k-1,r64) * Pvec(k-1)
          else
            Rvec(k) = -( 1.0_r64/(1.0_r64 - x) - 1.0_r64/(1.0_r64 + x) ) &
                      + real(k-1,r64) * Pvec(k-1)
          end if
        end if
      end do

      do j = 1, p
        do k = 1, p
          VT(k,j) = V(j,k)
        end do
      end do
      cr = Rvec
      call zge_r64(p, VT, cr)

      do j = 1, p
        Az_j    = cr(j) * (ci / (2.0_r64*pi*zsc))
        A1(j,m) = real(Az_j, r64)
        A2(j,m) = -aimag(Az_j)
      end do

      ! S = -(1/(1-x)^2 - (-1)^(k-1)/(1+x)^2)/2 + (k-1)*[0;R(1:p-1)]/2
      do k = 1, p
        if (k == 1) then
          Svec(k) = -( 1.0_r64/(1.0_r64 - x)**2 - 1.0_r64/(1.0_r64 + x)**2 ) / 2.0_r64
        else
          if (mod(k-1,2) == 0) then
            Svec(k) = -( 1.0_r64/(1.0_r64 - x)**2 - 1.0_r64/(1.0_r64 + x)**2 ) / 2.0_r64 &
                      + real(k-1,r64) * Rvec(k-1) / 2.0_r64
          else
            Svec(k) = -( 1.0_r64/(1.0_r64 - x)**2 + 1.0_r64/(1.0_r64 + x)**2 ) / 2.0_r64 &
                      + real(k-1,r64) * Rvec(k-1) / 2.0_r64
          end if
        end if
      end do

      do j = 1, p
        do k = 1, p
          VT(k,j) = V(j,k)
        end do
      end do
      cs2 = Svec
      call zge_r64(p, VT, cs2)

      do j = 1, p
        Azz_j   = cs2(j) * (ci / (2.0_r64*pi*zsc**2))
        A3(j,m) = real(Azz_j, r64)
        A4(j,m) = -aimag(Azz_j)
      end do

    end do

  end subroutine sdspecialquad_r64

  ! ------------------------------------------------------------------
  ! sdspecialquad_r128: same algorithm at quad precision.
  ! ------------------------------------------------------------------
  subroutine sdspecialquad_r128(nt, zt, p, zsrc, nzsrc, wzp, za, zb, iside, &
      As, Ad, A1, A2, A3, A4)
    integer(8),    intent(in)  :: nt, p
    complex(r128), intent(in)  :: zt(nt), zsrc(p), nzsrc(p), wzp(p), za, zb
    integer(8),    intent(in)  :: iside
    real(r128),    intent(out) :: As(p,nt), A1(p,nt), A2(p,nt), A3(p,nt), A4(p,nt)
    complex(r128), intent(out) :: Ad(p,nt)

    complex(r128) :: zsc, zmid, x, y(p), gam, log1mx, logprod
    complex(r128) :: c(p), V(p,p), VT(p,p)
    complex(r128) :: Pvec(p+1), Qvec(p), Rvec(p), Svec(p)
    complex(r128) :: cs(p), cd(p), cr(p), cs2(p)
    complex(r128) :: ci, one, Az_j, Azz_j
    real(r128)    :: pi, pi2inv, sign_k
    integer(8)    :: j, k, m

    pi     = acos(-1.0_r128)
    pi2inv = 1.0_r128 / (2.0_r128 * pi)
    ci     = cmplx(0.0_r128, 1.0_r128, r128)
    one    = cmplx(1.0_r128, 0.0_r128, r128)

    zsc  = (zb - za) / 2.0_r128
    zmid = (zb + za) / 2.0_r128

    do j = 1, p
      y(j) = (zsrc(j) - zmid) / zsc
    end do

    if (iside == 1) then
      gam = exp(-ci * pi / 4.0_r128)
    else
      gam = exp( ci * pi / 4.0_r128)
    end if

    do k = 1, p
      if (mod(k,2) == 1) then
        c(k) = cmplx(2.0_r128/real(k,r128), 0.0_r128, r128)
      else
        c(k) = cmplx(0.0_r128, 0.0_r128, r128)
      end if
    end do

    do j = 1, p
      V(j,1) = one
      do k = 2, p
        V(j,k) = V(j,k-1) * y(j)
      end do
    end do

    do m = 1, nt
      x = (zt(m) - zmid) / zsc

      log1mx = log(gam) + log((one - x) / (gam * (-one - x)))
      Pvec(1) = log1mx
      if (abs(x) <= 1.75_r128) then
        ! NEAR: upward recurrence
        do k = 1, p
          Pvec(k+1) = x * Pvec(k) + c(k)
        end do
      else
        ! FAR (|x|>1.75): seed top moment by direct quadrature, recurse DOWNWARD
        Pvec(p+1) = cmplx(0.0_r128, 0.0_r128, r128)
        do j = 1, p
          Pvec(p+1) = Pvec(p+1) + (wzp(j)/zsc) * (V(j,p)*y(j)) / (y(j) - x)
        end do
        do k = p, 2, -1
          Pvec(k) = (Pvec(k+1) - c(k)) / x
        end do
      end if

      logprod = log((one - x) * (-one - x))
      do k = 1, p
        if (mod(k,2) == 1) then
          Qvec(k) = (Pvec(k+1) - logprod) / real(k, r128)
        else
          Qvec(k) = (Pvec(k+1) - log1mx) / real(k, r128)
        end if
      end do

      do j = 1, p
        do k = 1, p
          VT(k,j) = V(j,k)
        end do
      end do
      cs = Qvec
      call zge_r128(p, VT, cs)
      do j = 1, p
        ! Match MATLAB SDspecialquad exactly: (1i*s.nx)' is conjugate-transpose.
        As(j,m) = real(cs(j) * conjg(ci*nzsrc(j)) * zsc, r128) * pi2inv &
                - log(abs(zsc)) * pi2inv * abs(wzp(j))
      end do

      do j = 1, p
        do k = 1, p
          VT(k,j) = V(j,k)
        end do
      end do
      cd = Pvec(1:p)
      call zge_r128(p, VT, cd)
      do j = 1, p
        Ad(j,m) = cd(j) * ci * pi2inv
      end do

      sign_k = 1.0_r128
      do k = 1, p
        if (k == 1) then
          Rvec(k) = -(one/(one-x) + sign_k*one/(one+x))
        else
          Rvec(k) = -(one/(one-x) + sign_k*one/(one+x)) &
                  + real(k-1, r128) * Pvec(k-1)
        end if
        sign_k = -sign_k
      end do
      do j = 1, p
        do k = 1, p
          VT(k,j) = V(j,k)
        end do
      end do
      cr = Rvec
      call zge_r128(p, VT, cr)
      do j = 1, p
        Az_j    = cr(j) * ci * pi2inv / zsc
        A1(j,m) =  real(Az_j, r128)
        A2(j,m) = -aimag(Az_j)
      end do

      sign_k = 1.0_r128
      do k = 1, p
        if (k == 1) then
          Svec(k) = (-one/(one-x)**2 + sign_k*one/(one+x)**2) / 2.0_r128
        else
          Svec(k) = (-one/(one-x)**2 + sign_k*one/(one+x)**2) / 2.0_r128 &
                  + real(k-1, r128) * Rvec(k-1) / 2.0_r128
        end if
        sign_k = -sign_k
      end do
      do j = 1, p
        do k = 1, p
          VT(k,j) = V(j,k)
        end do
      end do
      cs2 = Svec
      call zge_r128(p, VT, cs2)
      do j = 1, p
        Azz_j   = cs2(j) * ci * pi2inv / zsc**2
        A3(j,m) =  real(Azz_j, r128)
        A4(j,m) = -aimag(Azz_j)
      end do

    end do   ! m = 1..nt

  end subroutine sdspecialquad_r128

end module specialquad_mod
