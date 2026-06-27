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
  subroutine modal_green_r64(chi, n, vk, ve, Fn, An, dFn)
    real(r64),  intent(in)    :: chi
    integer(8), intent(in)    :: n
    real(r64),  intent(inout) :: vk, ve, Fn, An, dFn
    real(r64)  :: Kel, Eel, KmE, sm, Qm12, Q12, lam, vkF, QnpF, vkM, QnpM
    real(r64)  :: Qhi, Qlo, Qkm1, sc, sf, rk, twon1, mfac, mc, KmEmc, twopiinv
    real(r64)  :: Qf(0:n+1), F(0:n+1), dF(0:n+1)
    integer(8) :: k, M
    ! ---- vk = Q_{n-1/2}, ve ----
    sm = sqrt(2.0_r64/(chi+1.0_r64))
    call ellipke_mc_r64((chi-1.0_r64)/(chi+1.0_r64), Kel, Eel, KmE)
    Qm12 = sm*Kel; Q12 = chi*Qm12 - (chi+1.0_r64)*sm*Eel
    if (n == 0_8) then
      vk = Qm12; ve = (chi*Qm12 - Q12)/(chi+1.0_r64)
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
      lam = chi + sqrt(chi*chi - 1.0_r64)
      if (lam**(2.0_r64*real(n,r64)) < 1.0e6_r64) then
        vk = vkF; ve = twon1*(chi*vkF - QnpF)/(chi+1.0_r64)
      else
        vk = vkM; ve = twon1*(chi*vkM - QnpM)/(chi+1.0_r64)
      end if
    end if
    ! ---- Fn = P_{n-1/2}, An, dFn ----
    twopiinv = 2.0_r64/acos(-1.0_r64)
    mfac = 2.0_r64/(chi+1.0_r64); sm = sqrt(mfac); mc = (chi-1.0_r64)/(chi+1.0_r64)
    call ellipke_mc_r64(mfac, Kel, Eel, KmE)
    if (abs(mc) < 1.0e-300_r64) then
      KmEmc = acos(-1.0_r64)/4.0_r64
    else
      KmEmc = KmE/mc
    end if
    F(0) = sm*twopiinv*Kel
    F(1) = sm*twopiinv*((chi+1.0_r64)*Eel - Kel)
    do k = 1, n
      rk = real(k,r64)
      F(k+1) = (2.0_r64*rk*chi*F(k) - (rk-0.5_r64)*F(k-1))/(rk+0.5_r64)
    end do
    dF(0) = -sm*twopiinv*KmEmc/(2.0_r64*(chi+1.0_r64))
    dF(1) =  0.5_r64*sm*twopiinv*(-KmEmc/(chi+1.0_r64) + Eel)
    do k = 1, n
      rk = real(k,r64)
      dF(k+1) = (2.0_r64*rk*F(k) + 2.0_r64*rk*chi*dF(k) - (rk-0.5_r64)*dF(k-1))/(rk+0.5_r64)
    end do
    Fn = F(n); dFn = dF(n); An = -2.0_r64*(chi-1.0_r64)*dFn
  end subroutine modal_green_r64

  ! All-modes: fills vk(0:M)=Q_{n-1/2}, ve, Fn=P_{n-1/2}, An, dFn for n=0..M in ONE pass --
  ! single elliptic (AGM) eval per block and a single forward + single backward (Miller) recurrence sweep,
  ! replacing M+1 separate modal_green_r64 calls (each of which re-did the AGM and an O(6n) Miller sweep per mode).
  ! Same values as modal_green_r64 (low modes bit-identical via the forward recurrence; high modes machine-precision
  ! equal -- the Miller sweep here runs to Mtop(M) >= Mtop(n), i.e. at least as converged).
  subroutine modal_green_all_r64(chi, M, vk, ve, Fn, An, dFn)
    real(r64),  intent(in)    :: chi
    integer(8), intent(in)    :: M
    real(r64),  intent(inout) :: vk(0:M), ve(0:M), Fn(0:M), An(0:M), dFn(0:M)
    real(r64)  :: Kel, Eel, KmE, sm, Qm12, Q12, lam, sf, sc, rk, twon1, mfac, mc, KmEmc, twopiinv
    real(r64)  :: Qhi, Qlo, Qkm1
    real(r64)  :: Qf(0:M+1), Qm(0:M+1), F(0:M+1), dF(0:M+1)
    integer(8) :: k, n, Mtop
    ! ---- vk = Q_{n-1/2}, ve : forward recurrence (low modes) + ONE backward Miller sweep (high modes) ----
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
    ! ---- Fn = P_{n-1/2}, An, dFn : forward recurrence (dominant solution, stable) ----
    twopiinv = 2.0_r64/acos(-1.0_r64)
    mfac = 2.0_r64/(chi+1.0_r64); sm = sqrt(mfac); mc = (chi-1.0_r64)/(chi+1.0_r64)
    call ellipke_mc_r64(mfac, Kel, Eel, KmE)
    if (abs(mc) < 1.0e-300_r64) then
      KmEmc = acos(-1.0_r64)/4.0_r64
    else
      KmEmc = KmE/mc
    end if
    F(0) = sm*twopiinv*Kel
    F(1) = sm*twopiinv*((chi+1.0_r64)*Eel - Kel)
    do k = 1, M
      rk = real(k,r64)
      F(k+1) = (2.0_r64*rk*chi*F(k) - (rk-0.5_r64)*F(k-1))/(rk+0.5_r64)
    end do
    dF(0) = -sm*twopiinv*KmEmc/(2.0_r64*(chi+1.0_r64))
    dF(1) =  0.5_r64*sm*twopiinv*(-KmEmc/(chi+1.0_r64) + Eel)
    do k = 1, M
      rk = real(k,r64)
      dF(k+1) = (2.0_r64*rk*F(k) + 2.0_r64*rk*chi*dF(k) - (rk-0.5_r64)*dF(k-1))/(rk+0.5_r64)
    end do
    do n = 0, M
      Fn(n) = F(n); dFn(n) = dF(n); An(n) = -2.0_r64*(chi-1.0_r64)*dF(n)
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

end module axisym_modal_green_mod
