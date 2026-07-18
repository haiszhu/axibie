module axisym_modal_green_mod
  ! ------------------------------------------------------------------
  ! Axisymmetric modal Green's function carriers: the m-mode azimuthal Fourier
  ! component of the free-space Green's function and its pressure/derivative
  ! companions, for bodies of revolution.  In the static (kappa = 0,
  ! Laplace/Stokes) limit these are the half-integer-degree Legendre functions
  !   vk = Q_{n-1/2}(chi)            <- the modal Green's function itself,
  !   ve, Fn = P_{n-1/2}(chi), An, dFn   <- pressure / derivative companions.
  ! cf. Garritano, Kluger, Rokhlin & Serkh, J. Comput. Phys. 471 (2022) 111585
  ! ("efficient evaluation of the azimuthal Fourier components of the Green's
  ! function for Helmholtz's equation"), the general-kappa setting; this module
  ! is the kappa = 0 evaluator (forward-looking: a general-kappa routine can live
  ! here too).  Moved out of axistokes3d_mod (now only general numeric primitives);
  ! shared by both the Stokes and Laplace kernel-split + close-eval modules.
  ! ------------------------------------------------------------------
  use axistokes3d_mod, only: r64, ellipke_mc_r64
  implicit none

contains

  ! ------------------------------------------------------------------
  ! modal_green (single mode n)
  ! Half-integer Legendre carriers for the mode-n axisym close-eval split:
  !   vk = Q_{n-1/2}(chi),  ve = (2n+1)(chi Q_{n-1/2} - Q_{n+1/2})/(chi+1),
  !   Fn = P_{n-1/2}(chi),  An = (2n+1)(chi P_{n-1/2} - P_{n+1/2})/(chi+1),
  !   dFn = dP_{n-1/2}/dchi.
  ! Q-side (vk,ve): ellipke_mc at mc=(chi-1)/(chi+1), forward recurrence with
  ! scaled-Miller backward fallback (decaying-Q unstable for high n).  P-side
  ! (Fn,An,dFn): ellipke_mc at mc=2/(chi+1) (the complementary modulus); An,dFn
  ! are F-differences cancelling as chi->1, taken cancellation-free via K-E.
  ! ------------------------------------------------------------------
  subroutine modal_green_r64(chim1, n, vk, ve, Fn, An, dFn)
    real(r64),  intent(in)    :: chim1                 ! chi-1, carried to avoid the 1+small cancellation
    integer(8), intent(in)    :: n
    real(r64),  intent(inout) :: vk, ve, Fn, An, dFn
    real(r64)  :: chi, Kel, Eel, KmE, sm, Qm12, Q12, lam, vkF, QnpF, vkM, QnpM
    real(r64)  :: Qhi, Qlo, Qkm1, sc, sf, rk, twon1, mfac, mc, KmEmc, twopiinv
    real(r64)  :: Qf(0:n+1), F(0:n+1), dF(0:n+1)
    integer(8) :: k, M
    chi = 1.0_r64 + chim1                              ! reconstruct only for the (safe) recurrence multiplies
    ! ---- vk = Q_{n-1/2}, ve ----   mc=(chi-1)/(chi+1)=chim1/(chim1+2), chi+1=chim1+2 : all cancellation-free
    sm = sqrt(2.0_r64/(chim1+2.0_r64))
    call ellipke_mc_r64(chim1/(chim1+2.0_r64), Kel, Eel, KmE)
    Qm12 = sm*Kel; Q12 = chi*Qm12 - (chim1+2.0_r64)*sm*Eel
    if (n == 0_8) then
      vk = Qm12; ve = (chi*Qm12 - Q12)/(chim1+2.0_r64)
    else
      twon1 = 2.0_r64*real(n,r64) + 1.0_r64
      Qf(0) = Qm12; Qf(1) = Q12
      do k = 1, n
        rk = real(k,r64)
        Qf(k+1) = (2.0_r64*rk*chi*Qf(k) - (rk-0.5_r64)*Qf(k-1))/(rk+0.5_r64)
      end do
      vkF = Qf(n); QnpF = Qf(n+1)
      M = n + max(6_8*n, 80_8); sf = 2.0_r64**(-512)
      Qhi = 0.0_r64; Qlo = 1.0_r64; vkM = 0.0_r64; QnpM = 0.0_r64
      do k = M, 1, -1
        rk = real(k,r64)
        Qkm1 = (2.0_r64*rk*chi*Qlo - (rk+0.5_r64)*Qhi)/(rk-0.5_r64)
        if (k == n+1) QnpM = Qlo
        if (k == n)   vkM  = Qlo
        Qhi = Qlo; Qlo = Qkm1
        if (abs(Qlo) > 1.0e150_r64) then
          Qhi = Qhi*sf; Qlo = Qlo*sf; vkM = vkM*sf; QnpM = QnpM*sf
        end if
      end do
      sc = Qm12/Qlo; vkM = sc*vkM; QnpM = sc*QnpM
      lam = chi + sqrt(chim1*(chim1+2.0_r64))          ! chi^2-1 = chim1(chim1+2), cancellation-free
      if (lam**(2.0_r64*real(n,r64)) < 1.0e6_r64) then
        vk = vkF; ve = twon1*(chi*vkF - QnpF)/(chim1+2.0_r64)
      else
        vk = vkM; ve = twon1*(chi*vkM - QnpM)/(chim1+2.0_r64)
      end if
    end if
    ! ---- Fn = P_{n-1/2}, An, dFn ----
    twopiinv = 2.0_r64/acos(-1.0_r64)
    mfac = 2.0_r64/(chim1+2.0_r64); sm = sqrt(mfac); mc = chim1/(chim1+2.0_r64)
    call ellipke_mc_r64(mfac, Kel, Eel, KmE)
    if (abs(mc) < 1.0e-300_r64) then
      KmEmc = acos(-1.0_r64)/4.0_r64
    else
      KmEmc = KmE/mc
    end if
    F(0) = sm*twopiinv*Kel
    F(1) = sm*twopiinv*((chim1+2.0_r64)*Eel - Kel)
    do k = 1, n
      rk = real(k,r64)
      F(k+1) = (2.0_r64*rk*chi*F(k) - (rk-0.5_r64)*F(k-1))/(rk+0.5_r64)
    end do
    dF(0) = -sm*twopiinv*KmEmc/(2.0_r64*(chim1+2.0_r64))
    dF(1) =  0.5_r64*sm*twopiinv*(-KmEmc/(chim1+2.0_r64) + Eel)
    do k = 1, n
      rk = real(k,r64)
      dF(k+1) = (2.0_r64*rk*F(k) + 2.0_r64*rk*chi*dF(k) - (rk-0.5_r64)*dF(k-1))/(rk+0.5_r64)
    end do
    Fn = F(n); dFn = dF(n); An = -2.0_r64*chim1*dFn
  end subroutine modal_green_r64

  ! All-modes: fills vk(0:M)=Q_{n-1/2}, ve, Fn=P_{n-1/2}, An, dFn for n=0..M in ONE pass --
  ! single elliptic (AGM) eval per block and a single forward + single backward (Miller) recurrence sweep,
  ! replacing M+1 separate modal_green_r64 calls (each of which re-did the AGM and an O(6n) Miller sweep per mode).
  ! Same values as modal_green_r64 (low modes bit-identical via the forward recurrence; high modes machine-precision
  ! equal -- the Miller sweep here runs to Mtop(M) >= Mtop(n), i.e. at least as converged).
  subroutine modal_green_all_r64(chim1, M, vk, ve, Fn, An, dFn)
    real(r64),  intent(in)    :: chim1                 ! chi-1, carried to avoid the 1+small cancellation
    integer(8), intent(in)    :: M
    real(r64),  intent(inout) :: vk(0:M), ve(0:M), Fn(0:M), An(0:M), dFn(0:M)
    real(r64)  :: chi, Kel, Eel, KmE, sm, Qm12, Q12, lam, sf, sc, rk, twon1, mfac, mc, KmEmc, twopiinv
    real(r64)  :: Qhi, Qlo, Qkm1
    real(r64)  :: Qf(0:M+1), Qm(0:M+1), F(0:M+1), dF(0:M+1)
    integer(8) :: k, n, Mtop
    chi = 1.0_r64 + chim1                              ! reconstruct only for the (safe) recurrence multiplies
    ! ---- vk = Q_{n-1/2}, ve : forward recurrence (low modes) + ONE backward Miller sweep (high modes) ----
    sm = sqrt(2.0_r64/(chim1+2.0_r64))
    call ellipke_mc_r64(chim1/(chim1+2.0_r64), Kel, Eel, KmE)
    Qm12 = sm*Kel; Q12 = chi*Qm12 - (chim1+2.0_r64)*sm*Eel
    Qf(0) = Qm12; Qf(1) = Q12
    do k = 1, M
      rk = real(k,r64)
      Qf(k+1) = (2.0_r64*rk*chi*Qf(k) - (rk-0.5_r64)*Qf(k-1))/(rk+0.5_r64)
    end do
    Mtop = M + max(6_8*M, 80_8); sf = 2.0_r64**(-512)
    Qhi = 0.0_r64; Qlo = 1.0_r64; Qm = 0.0_r64
    do k = Mtop, 1, -1
      if (k <= M+1) Qm(k) = Qlo
      rk = real(k,r64)
      Qkm1 = (2.0_r64*rk*chi*Qlo - (rk+0.5_r64)*Qhi)/(rk-0.5_r64)
      Qhi = Qlo; Qlo = Qkm1
      if (abs(Qlo) > 1.0e150_r64) then
        Qhi = Qhi*sf; Qlo = Qlo*sf
        do n = max(k,0_8), M+1; Qm(n) = Qm(n)*sf; end do
      end if
    end do
    Qm(0) = Qlo; sc = Qm12/Qm(0)
    do k = 0, M+1; Qm(k) = sc*Qm(k); end do
    lam = chi + sqrt(chim1*(chim1+2.0_r64))            ! chi^2-1 = chim1(chim1+2), cancellation-free
    vk(0) = Qm12; ve(0) = (chi*Qm12 - Q12)/(chim1+2.0_r64)
    do n = 1, M
      twon1 = 2.0_r64*real(n,r64) + 1.0_r64
      if (lam**(2.0_r64*real(n,r64)) < 1.0e6_r64) then
        vk(n) = Qf(n); ve(n) = twon1*(chi*Qf(n) - Qf(n+1))/(chim1+2.0_r64)
      else
        vk(n) = Qm(n); ve(n) = twon1*(chi*Qm(n) - Qm(n+1))/(chim1+2.0_r64)
      end if
    end do
    ! ---- Fn = P_{n-1/2}, An, dFn : forward recurrence (dominant solution, stable) ----
    twopiinv = 2.0_r64/acos(-1.0_r64)
    mfac = 2.0_r64/(chim1+2.0_r64); sm = sqrt(mfac); mc = chim1/(chim1+2.0_r64)
    call ellipke_mc_r64(mfac, Kel, Eel, KmE)
    if (abs(mc) < 1.0e-300_r64) then
      KmEmc = acos(-1.0_r64)/4.0_r64
    else
      KmEmc = KmE/mc
    end if
    F(0) = sm*twopiinv*Kel
    F(1) = sm*twopiinv*((chim1+2.0_r64)*Eel - Kel)
    do k = 1, M
      rk = real(k,r64)
      F(k+1) = (2.0_r64*rk*chi*F(k) - (rk-0.5_r64)*F(k-1))/(rk+0.5_r64)
    end do
    dF(0) = -sm*twopiinv*KmEmc/(2.0_r64*(chim1+2.0_r64))
    dF(1) =  0.5_r64*sm*twopiinv*(-KmEmc/(chim1+2.0_r64) + Eel)
    do k = 1, M
      rk = real(k,r64)
      dF(k+1) = (2.0_r64*rk*F(k) + 2.0_r64*rk*chi*dF(k) - (rk-0.5_r64)*dF(k-1))/(rk+0.5_r64)
    end do
    do n = 0, M
      Fn(n) = F(n); dFn(n) = dF(n); An(n) = -2.0_r64*chim1*dF(n)
    end do
  end subroutine modal_green_all_r64

  ! All-modes FAR part: fills ONLY vk(0:M)=Q_{n-1/2} and ve(0:M) -- the Q-side first half --
  ! skipping the pressure/derivative companions Fn/An/dFn (and with them a second elliptic AGM
  ! eval and two forward recurrences).  The far-field block builders (SLP/DLP/SLPn/DLPn far
  ! codegen) use only vk, ve, so the P-side is wasted there.  Bit-identical vk/ve to
  ! modal_green_all_r64 (verbatim extraction of its first half).
  subroutine modal_green_all_far_r64(chi, M, vk, ve)
    real(r64),  intent(in)    :: chi
    integer(8), intent(in)    :: M
    real(r64),  intent(inout) :: vk(0:M), ve(0:M)
    real(r64)  :: Kel, Eel, KmE, sm, Qm12, Q12, lam, sf, sc, rk, twon1
    real(r64)  :: Qhi, Qlo, Qkm1
    real(r64)  :: Qf(0:M+1), Qm(0:M+1)
    integer(8) :: k, n, Mtop
    sm = sqrt(2.0_r64/(chi+1.0_r64))
    call ellipke_mc_r64((chi-1.0_r64)/(chi+1.0_r64), Kel, Eel, KmE)
    Qm12 = sm*Kel; Q12 = chi*Qm12 - (chi+1.0_r64)*sm*Eel
    Qf(0) = Qm12; Qf(1) = Q12
    do k = 1, M
      rk = real(k,r64)
      Qf(k+1) = (2.0_r64*rk*chi*Qf(k) - (rk-0.5_r64)*Qf(k-1))/(rk+0.5_r64)
    end do
    Mtop = M + max(6_8*M, 80_8); sf = 2.0_r64**(-512)
    Qhi = 0.0_r64; Qlo = 1.0_r64; Qm = 0.0_r64
    do k = Mtop, 1, -1
      if (k <= M+1) Qm(k) = Qlo
      rk = real(k,r64)
      Qkm1 = (2.0_r64*rk*chi*Qlo - (rk+0.5_r64)*Qhi)/(rk-0.5_r64)
      Qhi = Qlo; Qlo = Qkm1
      ! laplace green's identity test for target points near sphere north/south pole would reveal this
      if (abs(Qlo) > 1.0e150_r64) then
        Qhi = Qhi*sf; Qlo = Qlo*sf
        do n = max(k,0_8), M+1; Qm(n) = Qm(n)*sf; end do
      end if
    end do
    Qm(0) = Qlo; sc = Qm12/Qm(0)
    do k = 0, M+1; Qm(k) = sc*Qm(k); end do
    lam = chi + sqrt(chi*chi - 1.0_r64)
    vk(0) = Qm12; ve(0) = (chi*Qm12 - Q12)/(chi+1.0_r64)
    do n = 1, M
      twon1 = 2.0_r64*real(n,r64) + 1.0_r64
      if (lam**(2.0_r64*real(n,r64)) < 1.0e6_r64) then
        vk(n) = Qf(n); ve(n) = twon1*(chi*Qf(n) - Qf(n+1))/(chi+1.0_r64)
      else
        vk(n) = Qm(n); ve(n) = twon1*(chi*Qm(n) - Qm(n+1))/(chi+1.0_r64)
      end if
    end do
  end subroutine modal_green_all_far_r64

  ! Vectorized chunkie carrier: bit-for-bit port of chnk.axissymlap2d.qleg_half_miller_vec
  !   qm = Q_{k-1/2}(chi), qmd = Q'_{k-1/2}, qmdd = Q''_{k-1/2},  modes k=0..m, chi = t+1,
  ! one column per point, t(1:n) = chi-1 carried exactly as the MATLAB (never reconstructed
  ! from chi, so the near-surface log-seed sees the same bits).
  ! Mirrors the MATLAB structure: 3-regime seed [Q_{-1/2},Q_{1/2},Q'_{-1/2}] (near log-series /
  ! far 1/sqrt-series, element-wise; mid AGM run as ONE group with the shared max-residual
  ! break -- the only cross-point coupling in the MATLAB), then a per-point forward (chi<1.005)
  ! or backward-Miller recurrence (element-wise in the MATLAB, so per-point is bit-identical).
  ! Requires -ffp-contract=off (already in FFLAGS) for the exact chunkie match.
  subroutine modal_green_qleg_half_miller_vec_r64(n, t, m, qm, qmd, qmdd, pm, pmd, pmdd)
    integer(8), intent(in)    :: n, m
    real(r64),  intent(in)    :: t(n)                            ! chi - 1
    real(r64),  intent(inout) :: qm(0:m,n), qmd(0:m,n), qmdd(0:m,n)
    real(r64),  intent(inout), optional :: pm(0:m,n), pmd(0:m,n), pmdd(0:m,n)  ! P_{k-1/2},P',P'' (1st kind); absent -> Q-only
    real(r64), allocatable :: q0(:), q1(:)                       ! seeds Q_{-1/2}, Q_{1/2}, all points
    real(r64), allocatable :: agm_a(:), agm_a0(:), agm_b0(:), agm_aa0(:), agm_bb0(:), agm_fact(:)
    integer(8), allocatable :: imid(:)
    real(r64)  :: tj, chi, den, q0dj, b, v, x, tff, t0, pi
    real(r64)  :: delt, a1, b1, aa1, bb1, cc0, fF, fE, res, resmax
    real(r64)  :: fprev, fprevprev, f, dd, ratio
    real(r64)  :: mfac, mc, sm, twopiinv, Kel, Eel, KmE, KmEmc, p0, p1, p0d, p1d  ! P-side seed scalars
    real(r64)  :: cA(0:1002), cB(0:1002), cC(0:1002), cD(0:1002), cH(0:1002)  ! recurrence coef tables
    integer(8) :: i, j, jr, k, jj, nmid, it, nterms
    logical    :: use_fwd, fwd_allowed
    integer(8), parameter :: maxiter = 1000_8
    real(r64),  parameter :: upbound = 1.0e17_r64
    real(r64), parameter :: c0n(18) = (/ &                        ! qeval0_near (== qeval0der_near)
       1.7328679513998632735_r64, -0.5000000000000000000_r64, &
      -0.09160849392498290919_r64, 0.062500000000000000000_r64, &
       0.019905513916401443210_r64,-0.017578125000000000000_r64, &
      -0.006097834693194945559_r64, 0.006103515625000000000_r64, &
       0.0021674343381175963469_r64,-0.0023365020751953125000_r64, &
      -0.0008357538695841108955_r64, 0.0009462833404541015625_r64, &
       0.0003390851133264318725_r64,-0.00039757043123245239258_r64, &
      -0.00014242013770157621830_r64, 0.00017140153795480728149_r64, &
       0.00006133159221782055041_r64,-0.00007532294148404616863_r64 /)
    real(r64), parameter :: c1n(18) = (/ &                        ! qeval1_near
      -0.2671320486001367265_r64, -0.5000000000000000000_r64, &
       0.5248254817749487276_r64, -0.18750000000000000000_r64, &
      -0.04098835652733573868_r64, 0.029296875000000000000_r64, &
       0.009513531070472923783_r64,-0.008544921875000000000_r64, &
      -0.0029774361551467310174_r64, 0.0030040740966796875000_r64, &
       0.0010682069932178195667_r64,-0.0011565685272216796875_r64, &
      -0.0004138797762860294822_r64, 0.00046985596418380737305_r64, &
       0.00016838776925222835322_r64,-0.00019777100533246994019_r64, &
      -0.00007084821236213522235_r64, 0.00008536600034858565778_r64 /)
    real(r64), parameter :: c0f(9) = (/ &                         ! qeval0_far
       2.2214414690791831235_r64, 0.0_r64, 0.41652027545234683566_r64, 0.0_r64, &
       0.22778452563800217575_r64, 0.0_r64, 0.15660186137612649583_r64, 0.0_r64, &
       0.11928657409509635424_r64 /)
    real(r64), parameter :: c1f(9) = (/ &                         ! qeval1_far
       0.55536036726979578088_r64, 0.0_r64, 0.26032517215771677229_r64, 0.0_r64, &
       0.17083839422850163181_r64, 0.0_r64, 0.12723901236810277786_r64, 0.0_r64, &
       0.10139358798083190111_r64 /)

    pi = acos(-1.0_r64)

    ! ---- recurrence coefficient tables, shared by every point (the same single divisions the
    !      vectorized MATLAB evaluates once per iteration; downstream products bit-identical) ----
    cA(0:1) = 0.0_r64; cB(0:1) = 0.0_r64; cC(0:1) = 0.0_r64; cD(0:1) = 0.0_r64
    do jr = 2, 1002
      cA(jr) = 4.0_r64*real(jr-1,r64)/real(2*jr-1,r64)   ! forward form   4(j-1)/(2j-1)
      cB(jr) = real(2*jr-3,r64)/real(2*jr-1,r64)         !                (2j-3)/(2j-1)
      cC(jr) = 4.0_r64*real(jr-1,r64)/real(2*jr-3,r64)   ! backward form  4(j-1)/(2j-3)
      cD(jr) = real(2*jr-1,r64)/real(2*jr-3,r64)         !                (2j-1)/(2j-3)
    end do
    do k = 0, 1002
      cH(k) = real(k,r64) - 0.5_r64                      ! k - 1/2
    end do
    fwd_allowed = (m <= 163_8)                           ! axmg thresholds are always-false for m>163

    ! ---- seed, mid regime (0.01 <= t <= 100) FIRST: ONE AGM over the whole mid subset with
    !      the shared max-residual break (the MATLAB's only cross-point coupling) ----
    allocate(q0(n), q1(n), imid(n))
    nmid = 0
    do j = 1, n
      if (t(j) >= 0.01_r64 .and. t(j) <= 100.0_r64) then
        nmid = nmid + 1; imid(nmid) = j
      end if
    end do
    if (nmid > 0) then
      allocate(agm_a(nmid), agm_a0(nmid), agm_b0(nmid), agm_aa0(nmid), agm_bb0(nmid), agm_fact(nmid))
      do jj = 1, nmid
        x = 1.0_r64 + t(imid(jj))
        agm_a(jj)   = sqrt(2.0_r64/(x+1.0_r64))
        delt        = 1.0_r64/sqrt(1.0_r64 - agm_a(jj)*agm_a(jj))
        agm_aa0(jj) = delt + sqrt(delt*delt - 1.0_r64)
        agm_bb0(jj) = 1.0_r64/(delt + sqrt(delt*delt - 1.0_r64))
        agm_a0(jj)  = 1.0_r64
        agm_b0(jj)  = 1.0_r64/delt
        agm_fact(jj) = ((agm_a0(jj)+agm_b0(jj))/2.0_r64)**2
      end do
      do it = 1, 1000
        resmax = 0.0_r64
        do jj = 1, nmid
          a1  = (agm_a0(jj)+agm_b0(jj))/2.0_r64;   b1  = sqrt(agm_a0(jj)*agm_b0(jj))
          aa1 = (agm_aa0(jj)+agm_bb0(jj))/2.0_r64; bb1 = sqrt(agm_aa0(jj)*agm_bb0(jj))
          agm_a0(jj) = a1; agm_b0(jj) = b1; agm_aa0(jj) = aa1; agm_bb0(jj) = bb1
          cc0 = (a1-b1)/2.0_r64
          agm_fact(jj) = agm_fact(jj) - cc0*cc0*2.0_r64**it
          res = abs(agm_a0(jj)-agm_b0(jj))/abs(agm_a0(jj)) + abs(agm_aa0(jj)-agm_bb0(jj))/abs(agm_aa0(jj))
          if (res > resmax) resmax = res
        end do
        if (resmax < 2.0_r64*1.0e-15_r64) exit
      end do
      do jj = 1, nmid
        j = imid(jj); x = 1.0_r64 + t(j)
        fF = pi/(2.0_r64*agm_aa0(jj)*sqrt(1.0_r64 - agm_a(jj)*agm_a(jj)))
        fE = pi*agm_fact(jj)/(2.0_r64*agm_a0(jj))
        q0(j) = agm_a(jj)*fF
        q1(j) = x*q0(j) - 2.0_r64*fE/agm_a(jj)
      end do
      deallocate(agm_a, agm_a0, agm_b0, agm_aa0, agm_bb0, agm_fact)
    end if
    deallocate(imid)

    ! ---- everything else is element-wise in the MATLAB -> one per-point pass, points are
    !      independent so the loop is OpenMP-parallel (dynamic schedule: near-cutoff Miller
    !      points run ~25x longer than far ones); results are bit-identical to serial ----
    !!$omp parallel do default(shared) schedule(dynamic,64) &
    !!$omp&  private(j, tj, chi, den, q0dj, b, v, x, tff, t0, i, jr, k, &
    !!$omp&          fprev, fprevprev, f, dd, ratio, nterms, use_fwd)
    do j = 1, n
      tj  = t(j)
      chi = tj + 1.0_r64                            ! MATLAB: chi = t+1
      den = -tj*tj - 2.0_r64*tj                     ! recurrence denominator -(t^2+2t)

      if (tj < 0.01_r64) then                       ! near: log-series (incl. qeval0der_near)
        b = log(tj)
        q0(j) = 0.0_r64; q1(j) = 0.0_r64; v = 1.0_r64
        do i = 1, 9
          q0(j) = q0(j) + v*c0n(2*i-1) + v*b*c0n(2*i)
          q1(j) = q1(j) + v*c1n(2*i-1) + v*b*c1n(2*i)
          v  = v*tj
        end do
        q0dj = c0n(2)/tj; v = 1.0_r64
        do i = 2, 9
          q0dj = q0dj + (i-1)*(v*c0n(2*i-1) + v*b*c0n(2*i)) + v*c0n(2*i)
          v    = v*tj
        end do
      else if (tj > 100.0_r64) then                 ! far: 1/sqrt series
        x = 1.0_r64 + tj; tff = 1.0_r64/x
        q0(j) = 0.0_r64; t0 = 1.0_r64/sqrt(x)
        do i = 1, 9
          q0(j) = q0(j) + t0*c0f(i); t0 = t0*tff
        end do
        q1(j) = 0.0_r64; t0 = 1.0_r64/sqrt(x)**3
        do i = 1, 9
          q1(j) = q1(j) + t0*c1f(i); t0 = t0*tff
        end do
        q0dj = (-0.5_r64*q1(j) + 0.5_r64*(1.0_r64+tj)*q0(j))/(-tj*(tj+2.0_r64))   ! seed: -t(t+2)
      else                                          ! mid: q0,q1 from the group AGM above
        q0dj = (-0.5_r64*q1(j) + 0.5_r64*(1.0_r64+tj)*q0(j))/(-tj*(tj+2.0_r64))
      end if

      use_fwd = fwd_allowed .and. (chi < 1.005_r64)

      if (use_fwd) then                             ! ===== FORWARD =====
        qm(0,j)  = q0(j); qm(1,j) = q1(j)
        qmd(0,j) = q0dj
        qmd(1,j) = (-q0(j) + chi*q1(j))/(2.0_r64*(chi+1.0_r64)*tj)
        do i = 1, m-1
          jr = i+1
          qm(i+1,j)  = cA(jr)*chi*qm(i,j) - cB(jr)*qm(i-1,j)
          qmd(i+1,j) = (cH(jr)*qm(i,j) - cH(jr)*chi*qm(i+1,j))/den
        end do
      else                                          ! ===== BACKWARD (Miller) =====
        fprev = 1.0_r64; fprevprev = 0.0_r64        ! 1. forward until it blows up -> nterms
        nterms = 0_8
        do i = m, maxiter
          jr = i+1
          f  = cA(jr)*chi*fprev - cB(jr)*fprevprev
          dd = (cH(jr)*fprevprev - cH(jr)*chi*fprev)/den
          if (abs(f) >= upbound .and. abs(dd) >= upbound) then
            nterms = i+1; exit
          end if
          fprevprev = fprev; fprev = f
        end do
        if (nterms == 0_8) nterms = maxiter
        fprev = 1.0_r64; fprevprev = 0.0_r64        ! 2. backward from nterms to m
        do i = 1, nterms - m + 1
          jr = nterms - i + 2
          f = cC(jr)*chi*fprev - cD(jr)*fprevprev
          if (jr <= nterms) then
            fprevprev = fprev; fprev = f
          else
            fprevprev = 0.0_r64
          end if
        end do
        qm(m-1,j) = fprev; qm(m,j) = fprevprev
        do i = 1, m-1                               ! backward from m to 1
          jr = m - i + 1
          qm(jr-2,j) = cC(jr)*chi*qm(jr-1,j) - cD(jr)*qm(jr,j)
          ! laplace green's identity test for target points near sphere north/south pole would reveal this
          if (abs(qm(jr-2,j)) > 1.0e150_r64) then
            do k = jr-2, m; qm(k,j) = qm(k,j)*2.0_r64**(-512); end do
          end if
        end do
        ratio = q0(j) / qm(0,j)                     ! 3. scale by known Q_{-1/2}
        do k = 0, m
          qm(k,j) = qm(k,j)*ratio
        end do
        qmd(0,j) = q0dj                             ! 4. first derivatives (backward form)
        do i = 1, m
          qmd(i,j) = (-cH(i)*chi*qm(i,j) + cH(i)*qm(i-1,j))/den
        end do
      end if

      do k = 0, m                                   ! 5. second derivatives (Legendre ODE)
        qmdd(k,j) = (-cH(k)*cH(k+1)*qm(k,j) + 2.0_r64*chi*qmd(k,j))/den
      end do
    end do
    !!$omp end parallel do

    ! ================= P-side: P_{k-1/2}, P', P'' (1st kind, DOMINANT -> always FORWARD) =================
    !   OPTIONAL outputs (pm/pmd/pmdd absent -> Q-only, the Q block above is untouched).  P solves the SAME
    !   recurrence/ODE as Q but is the dominant solution, so plain forward is unconditionally stable (no
    !   Miller): seed [P_{-1/2},P_{1/2},P'_{-1/2},P'_{1/2}] from K(mc),E(mc) at the complementary modulus
    !   mc = t/(t+2) via ellipke_mc_r64(mfac), mfac = 2/(chi+1) (accurate as chi->1, where mc->0), then the
    !   SAME cA/cB forward + cH derivative-identity + Legendre-ODE recurrences as the Q forward block.
    !   Mirrors the P-side of axmg_modal_green_qleg_half_miller_vec.m.
    if (present(pm)) then
      twopiinv = 2.0_r64/pi
      !!$omp parallel do default(shared) schedule(dynamic,64) &
      !!$omp&  private(j, tj, chi, den, mfac, mc, sm, Kel, Eel, KmE, KmEmc, p0, p1, p0d, p1d, i, jr, k)
      do j = 1, n
        tj  = t(j)
        chi = tj + 1.0_r64                           ! MATLAB: chip = t+1
        den = -tj*tj - 2.0_r64*tj                    ! -(t^2+2t)
        mfac = 2.0_r64/(tj+2.0_r64)                  ! m  = 2/(chi+1)
        mc   = tj/(tj+2.0_r64)                       ! mc = (chi-1)/(chi+1)
        sm   = sqrt(mfac)
        call ellipke_mc_r64(mfac, Kel, Eel, KmE)     ! K(mc), E(mc), K-E   (mc = 1 - mfac)
        if (mc < 1.0e-300_r64) then
          KmEmc = pi/4.0_r64
        else
          KmEmc = KmE/mc
        end if
        p0  = sm*twopiinv*Kel                                     ! P_{-1/2}
        p1  = sm*twopiinv*((tj+2.0_r64)*Eel - Kel)               ! P_{1/2}
        p0d = -sm*twopiinv*KmEmc/(2.0_r64*(tj+2.0_r64))          ! P'_{-1/2}
        p1d = (-p0 + chi*p1)/(2.0_r64*(tj+2.0_r64)*tj)           ! P'_{1/2}  (same identity as q1d)
        pm(0,j)  = p0;  pm(1,j)  = p1
        pmd(0,j) = p0d; pmd(1,j) = p1d
        do i = 1, m-1                                ! forward (dominant -> stable, no Miller)
          jr = i+1
          pm(i+1,j)  = cA(jr)*chi*pm(i,j) - cB(jr)*pm(i-1,j)
          pmd(i+1,j) = (cH(jr)*pm(i,j) - cH(jr)*chi*pm(i+1,j))/den
        end do
        do k = 0, m                                  ! second derivatives (Legendre ODE)
          pmdd(k,j) = (-cH(k)*cH(k+1)*pm(k,j) + 2.0_r64*chi*pmd(k,j))/den
        end do
      end do
      !!$omp end parallel do
    end if
  end subroutine modal_green_qleg_half_miller_vec_r64

end module axisym_modal_green_mod
