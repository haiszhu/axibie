module axissymstok_specialquad_mod
  use axistokes3d_mod, only: r64, gauss_r64, carrier_r64, lagrange_interp_r64
  use specialquad_mod, only: sdspecialquad_r64
  use axissymstok_kernelsplit_mod, only: axissymstok_slp_coef_r64, axissymstok_slp_coef_nmode_r64, &
       axissymstok_slpn_coef_r64, axissymstok_slpn_coef_nmode_r64, &
       axissymstok_dlp_coef_r64, axissymstok_dlp_coef_nmode_r64, weight_setup_r64, &
       axissymstok_dlp_aziquad_r64, axissymstok_slpn_aziquad_r64, axissymstok_dlpn_coef_r64, &
       axissymstok_dlpn_coef_nmode_r64, axissymstok_dlpn_aziquad_r64
  implicit none
  private
  public :: axissymstok_slp_blockmat_r64
  public :: axissymstok_slp_blockmat_nmode_r64
  public :: axissymstok_slpn_blockmat_r64
  public :: axissymstok_slpn_blockmat_nmode_r64
  public :: axissymstok_dlp_blockmat_r64
  public :: axissymstok_dlp_blockmat_nmode_r64
  public :: axissymstok_dlpn_blockmat_r64
  public :: axissymstok_dlpn_blockmat_nmode_r64
contains

  subroutine axissymstok_slp_blockmat_r64(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, mu, A)
    integer(8),   intent(in)    :: nt, p, np, iside, iclosed
    complex(r64), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: sws(p*np), tpan(np+1), mu
    real(r64),    intent(inout) :: A(2*nt, 2*np*p)
    integer(8) :: q, nso, i, j, jp, iq, ia, l, nc, cjo, k, jo, ne, npa
    real(r64)  :: muinv, twopi, ipi, sumws, rt, zti, rho, rhop, zh, drho, rr2, chi, pref, spd
    real(r64)  :: SKrr, SErr, SKrz, SErz, SKzr, SEzr, SKzz, SEzz, vk, ve, Fn, An, dFn
    real(r64)  :: t1, tN, tm, denom, rlo, rhi, tgi, split
    complex(r64) :: ic, za, zb
    real(r64)    :: tglp(p), wglp(p), Dp(p,p), tglq(2*p), wglq(2*p), Dq(2*p,2*p), IP2(2*p,p)
    real(r64)    :: tc(p,np), rk(p), re2(2), IPk(p,p), IPe2(2,p), wsq(2*p), wsp(p)
    complex(r64) :: Yp(p), Yq(2*p), dYp(p), dYq(2*p), nvq(2*p), wxpq(2*p), ec2(2), xc(p,np)
    real(r64),    allocatable :: tin(:)
    logical,      allocatable :: cl(:)
    integer(8),   allocatable :: ci(:)
    complex(r64), allocatable :: zc(:)
    real(r64),    allocatable :: C1(:,:), C2(:,:), C3(:,:), As(:,:), A1(:,:), A2(:,:), A3(:,:), A4(:,:)
    real(r64),    allocatable :: Gcl(:,:), Gp(:,:), ff(:,:)
    complex(r64), allocatable :: Ad(:,:)
    q = 2*p; nso = np*p; ic = (0.0_r64,1.0_r64); twopi = 2.0_r64*acos(-1.0_r64)
    ipi = 1.0_r64/acos(-1.0_r64); muinv = 1.0_r64/mu
    call gauss_r64(p,     tglp, wglp, Dp)
    call gauss_r64(2_8*p, tglq, wglq, Dq)
    call lagrange_interp_r64(p, tglp, q, tglq, IP2)
    do j = 1, np
      do i = 1, p
        xc(i,j) = sx((j-1)*p+i)
        tc(i,j) = tpan(j) + (1.0_r64+tglp(i))/2.0_r64*(tpan(j+1)-tpan(j))
      end do
    end do
    t1 = tc(1,1); tN = tc(p,np)
    ! ---- dyadic pole refinement of tpan -> tin (mirror StoSLPAxiBlockMat_v0) ----
    allocate(tin(np+1+256)); ne = np+1; tin(1:ne) = tpan(1:np+1)
    do while (tin(2) > t1)
      split = 0.5_r64*(tin(1)+tin(2))
      do i = ne, 2, -1; tin(i+1) = tin(i); end do
      tin(2) = split; ne = ne+1
    end do
    do while (tin(ne-1) < tN)
      split = 0.5_r64*(tin(ne)+tin(ne-1))
      tin(ne+1) = tin(ne); tin(ne) = split; ne = ne+1
    end do
    npa = ne-1
    allocate(cl(nt), ci(nt), zc(nt), C1(2*nt,2*q), C2(2*nt,2*q), C3(2*nt,2*q))
    allocate(As(q,nt), Ad(q,nt), A1(q,nt), A2(q,nt), A3(q,nt), A4(q,nt), Gcl(2*nt,2*q), Gp(nt,p), ff(p,4))
    A = 0.0_r64
    do k = 1, npa
      ! ---- this refined panel: parent jo, coarse->refined interp IPk, positions/endpoints ----
      tm = 0.5_r64*(tin(k)+tin(k+1)); jo = 1
      do j = 1, np
        if (tm >= tpan(j) .and. tm < tpan(j+1)) jo = j
      end do
      denom = tpan(jo+1) - tpan(jo)
      do i = 1, p
        tgi = tin(k) + (1.0_r64+tglp(i))/2.0_r64*(tin(k+1)-tin(k))
        rk(i) = (2.0_r64*tgi - (tpan(jo)+tpan(jo+1)))/denom
      end do
      call lagrange_interp_r64(p, tglp, p, rk, IPk)
      Yp = matmul(IPk, xc(:,jo))
      rlo = (2.0_r64*tin(k)   - (tpan(jo)+tpan(jo+1)))/denom
      rhi = (2.0_r64*tin(k+1) - (tpan(jo)+tpan(jo+1)))/denom
      re2(1) = rlo; re2(2) = rhi
      call lagrange_interp_r64(p, tglp, 2_8, re2, IPe2)
      ec2 = matmul(IPe2, xc(:,jo)); za = ec2(1); zb = ec2(2)
      cjo = (jo-1)*p
      ! ---- refined-panel geometry (upsample + speed/normal via D-matrix) ----
      Yq = matmul(IP2, Yp); dYq = matmul(Dq, Yq)
      do iq = 1, q
        spd = abs(dYq(iq)); nvq(iq) = -ic*dYq(iq)/spd; wsq(iq) = wglq(iq)*spd; wxpq(iq) = dYq(iq)*wglq(iq)
      end do
      dYp = matmul(Dp, Yp)
      do jp = 1, p
        wsp(jp) = wglp(jp)*abs(dYp(jp))
      end do
      sumws = sum(wsp)
      nc = 0
      do i = 1, nt
        cl(i) = (abs(tx(i)-za) + abs(tx(i)-zb)) < 1.5_r64*sumws
        if (cl(i)) then; nc = nc + 1; ci(nc) = i; zc(nc) = tx(i); end if
      end do
      ! ---- near: 2p close-eval (coef + sdspecialquad), fold q->p (IP2) -> coarse (IPk) ----
      if (nc > 0) then
        call axissymstok_slp_coef_r64(nc, zc(1:nc), q, Yq, mu, C1(1:2*nc,1:2*q), C2(1:2*nc,1:2*q), C3(1:2*nc,1:2*q))
        call sdspecialquad_r64(nc, zc(1:nc), q, Yq, nvq, wxpq, za, zb, iside, &
             As(1:q,1:nc), Ad(1:q,1:nc), A1(1:q,1:nc), A2(1:q,1:nc), A3(1:q,1:nc), A4(1:q,1:nc))
        do iq = 1, q
          do ia = 1, nc
            Gcl(ia,      iq)     = twopi*C1(ia,iq)*As(iq,ia)         + C2(ia,iq)*wsq(iq) &
                                 + twopi*real(Ad(iq,ia)*(C3(ia,iq)/nvq(iq)), r64)
            Gcl(ia,      q+iq)   = twopi*C1(ia,q+iq)*As(iq,ia)       + C2(ia,q+iq)*wsq(iq) &
                                 + twopi*real(Ad(iq,ia)*(C3(ia,q+iq)/nvq(iq)), r64)
            Gcl(nc+ia,   iq)     = twopi*C1(nc+ia,iq)*As(iq,ia)      + C2(nc+ia,iq)*wsq(iq) &
                                 + twopi*real(Ad(iq,ia)*(C3(nc+ia,iq)/nvq(iq)), r64)
            Gcl(nc+ia,   q+iq)   = twopi*C1(nc+ia,q+iq)*As(iq,ia)    + C2(nc+ia,q+iq)*wsq(iq) &
                                 + twopi*real(Ad(iq,ia)*(C3(nc+ia,q+iq)/nvq(iq)), r64)
          end do
        end do
        Gp(1:nc,:) = matmul(Gcl(1:nc, 1:q), IP2)
        do ia = 1, nc; do l = 1, p; A(ci(ia),    cjo+l)     = A(ci(ia),    cjo+l)     + sum(Gp(ia,1:p)*IPk(:,l)); end do; end do
        Gp(1:nc,:) = matmul(Gcl(1:nc, q+1:2*q), IP2)
        do ia = 1, nc; do l = 1, p; A(ci(ia),    nso+cjo+l) = A(ci(ia),    nso+cjo+l) + sum(Gp(ia,1:p)*IPk(:,l)); end do; end do
        Gp(1:nc,:) = matmul(Gcl(nc+1:2*nc, 1:q), IP2)
        do ia = 1, nc; do l = 1, p; A(nt+ci(ia), cjo+l)     = A(nt+ci(ia), cjo+l)     + sum(Gp(ia,1:p)*IPk(:,l)); end do; end do
        Gp(1:nc,:) = matmul(Gcl(nc+1:2*nc, q+1:2*q), IP2)
        do ia = 1, nc; do l = 1, p; A(nt+ci(ia), nso+cjo+l) = A(nt+ci(ia), nso+cjo+l) + sum(Gp(ia,1:p)*IPk(:,l)); end do; end do
      end if
      ! ---- far: analytic 0th SLP (K,E two-carrier) on the refined p nodes, fold via IPk ----
      do i = 1, nt
        if (cl(i)) cycle
        rt = real(tx(i),r64); zti = aimag(tx(i))
        do jp = 1, p
          rhop = real(Yp(jp),r64); zh = zti - aimag(Yp(jp)); rho = rt; drho = rho - rhop
          rr2 = drho*drho + zh*zh; chi = 1.0_r64 + rr2/(2.0_r64*rho*rhop)
          call carrier_r64(chi, 0_8, vk, ve, Fn, An, dFn)
          pref = sqrt(rhop/rho)*ipi/8.0_r64 * wsp(jp) * muinv
          SKrr = (rho*rho + rhop*rhop + 2.0_r64*zh*zh)/(rho*rhop); SErr = -2.0_r64-4.0_r64*chi+2.0_r64*chi*drho*drho/rr2
          SKrz = zh/rho;   SErz = 2.0_r64*(drho*zh/rr2 - zh/(2.0_r64*rho))
          SKzr = -zh/rhop; SEzr = 2.0_r64*(drho*zh/rr2 + zh/(2.0_r64*rhop))
          SKzz = 2.0_r64;  SEzz = 2.0_r64*zh*zh/rr2
          ff(jp,1) = pref*(SKrr*vk + SErr*ve); ff(jp,2) = pref*(SKrz*vk + SErz*ve)
          ff(jp,3) = pref*(SKzr*vk + SEzr*ve); ff(jp,4) = pref*(SKzz*vk + SEzz*ve)
        end do
        do l = 1, p
          A(i,     cjo+l)     = A(i,     cjo+l)     + sum(ff(1:p,1)*IPk(:,l))
          A(i,     nso+cjo+l) = A(i,     nso+cjo+l) + sum(ff(1:p,2)*IPk(:,l))
          A(nt+i,  cjo+l)     = A(nt+i,  cjo+l)     + sum(ff(1:p,3)*IPk(:,l))
          A(nt+i,  nso+cjo+l) = A(nt+i,  nso+cjo+l) + sum(ff(1:p,4)*IPk(:,l))
        end do
      end do
    end do
    deallocate(tin, cl, ci, zc, C1, C2, C3, As, Ad, A1, A2, A3, A4, Gcl, Gp, ff)
  end subroutine axissymstok_slp_blockmat_r64

  subroutine axissymstok_slp_blockmat_nmode_r64(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, mu, A)
    integer(8),   intent(in)    :: nt, p, np, M, iside, iclosed
    complex(r64), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: sws(p*np), tpan(np+1), mu
    complex(r64), intent(inout) :: A(3*nt, 3*np*p, M+1)
    integer(8) :: q, nso, i, j, jp, iq, ia, l, nc, cjo, e, md, ar, b, k, jo, ne, npa
    real(r64)  :: ipi, twopi, muinv, sumws, rho, rhop, zh, rr2, chi, ws, spd, rt, zti, vk, ve, Fn, An, dFn
    real(r64)  :: rn, n2, fn2, cm1, cp1, rr, rrt, srt, t1, tN, tm, denom, rlo, rhi, tgi, split
    complex(r64) :: ic, za, zb, cc1, cc2, cc3
    complex(r64) :: sk1,se1,sk2,se2,sk3,se3,sk4,se4,sk5,se5,sk6,sk7,se7,sk8,sk9,se9
    real(r64)    :: tglp(p), wglp(p), Dp(p,p), tglq(2*p), wglq(2*p), Dq(2*p,2*p), IP2(2*p,p)
    real(r64)    :: tc(p,np), rk(p), re2(2), IPk(p,p), IPe2(2,p), wsq(2*p), wsp(p)
    complex(r64) :: Yp(p), Yq(2*p), dYp(p), dYq(2*p), nvq(2*p), wxpq(2*p), ec2(2), xc(p,np)
    real(r64),    allocatable :: tin(:)
    logical,      allocatable :: cl(:)
    integer(8),   allocatable :: ci(:)
    complex(r64), allocatable :: zc(:), C1(:,:,:), C2(:,:,:), C3(:,:,:), Gcq(:,:), Gpc(:,:), ff(:,:)
    real(r64),    allocatable :: As(:,:), A1(:,:), A2(:,:), A3(:,:), A4(:,:)
    complex(r64), allocatable :: Ad(:,:)
    q = 2*p; nso = np*p; ic = (0.0_r64,1.0_r64); ipi = 1.0_r64/acos(-1.0_r64)
    twopi = 2.0_r64*acos(-1.0_r64); muinv = 1.0_r64/mu
    call gauss_r64(p,     tglp, wglp, Dp)
    call gauss_r64(2_8*p, tglq, wglq, Dq)
    call lagrange_interp_r64(p, tglp, q, tglq, IP2)
    do j = 1, np
      do i = 1, p
        xc(i,j) = sx((j-1)*p+i)
        tc(i,j) = tpan(j) + (1.0_r64+tglp(i))/2.0_r64*(tpan(j+1)-tpan(j))
      end do
    end do
    t1 = tc(1,1); tN = tc(p,np)
    ! ---- dyadic pole refinement of tpan -> tin (mirror StoSLPAxiBlockMat_nmode) ----
    allocate(tin(np+1+256)); ne = np+1; tin(1:ne) = tpan(1:np+1)
    do while (tin(2) > t1)
      split = 0.5_r64*(tin(1)+tin(2))
      do i = ne, 2, -1; tin(i+1) = tin(i); end do
      tin(2) = split; ne = ne+1
    end do
    do while (tin(ne-1) < tN)
      split = 0.5_r64*(tin(ne)+tin(ne-1))
      tin(ne+1) = tin(ne); tin(ne) = split; ne = ne+1
    end do
    npa = ne-1
    allocate(cl(nt), ci(nt), zc(nt), C1(3*nt,3*q,M+1), C2(3*nt,3*q,M+1), C3(3*nt,3*q,M+1))
    allocate(As(q,nt), Ad(q,nt), A1(q,nt), A2(q,nt), A3(q,nt), A4(q,nt), Gcq(nt,q), Gpc(nt,p), ff(p,9))
    A = (0.0_r64,0.0_r64)
    do k = 1, npa
      ! ---- this refined panel: parent jo, coarse->refined interp IPk, positions/endpoints ----
      tm = 0.5_r64*(tin(k)+tin(k+1)); jo = 1
      do j = 1, np
        if (tm >= tpan(j) .and. tm < tpan(j+1)) jo = j
      end do
      denom = tpan(jo+1) - tpan(jo)
      do i = 1, p
        tgi = tin(k) + (1.0_r64+tglp(i))/2.0_r64*(tin(k+1)-tin(k))
        rk(i) = (2.0_r64*tgi - (tpan(jo)+tpan(jo+1)))/denom
      end do
      call lagrange_interp_r64(p, tglp, p, rk, IPk)
      Yp = matmul(IPk, xc(:,jo))
      rlo = (2.0_r64*tin(k)   - (tpan(jo)+tpan(jo+1)))/denom
      rhi = (2.0_r64*tin(k+1) - (tpan(jo)+tpan(jo+1)))/denom
      re2(1) = rlo; re2(2) = rhi
      call lagrange_interp_r64(p, tglp, 2_8, re2, IPe2)
      ec2 = matmul(IPe2, xc(:,jo)); za = ec2(1); zb = ec2(2)
      cjo = (jo-1)*p
      ! ---- refined-panel geometry (upsample + speed/normal via D-matrix) ----
      Yq = matmul(IP2, Yp); dYq = matmul(Dq, Yq)
      do iq = 1, q
        spd = abs(dYq(iq)); nvq(iq) = -ic*dYq(iq)/spd; wsq(iq) = wglq(iq)*spd; wxpq(iq) = dYq(iq)*wglq(iq)
      end do
      dYp = matmul(Dp, Yp)
      do jp = 1, p
        wsp(jp) = wglp(jp)*abs(dYp(jp))
      end do
      sumws = sum(wsp)
      nc = 0
      do i = 1, nt
        cl(i) = (abs(tx(i)-za) + abs(tx(i)-zb)) < 1.5_r64*sumws
        if (cl(i)) then; nc = nc + 1; ci(nc) = i; zc(nc) = tx(i); end if
      end do
      ! ---- near: all-mode SLP close (coef + sdspecialquad), fold q->p_coarse via IPqc ----
      if (nc > 0) then
        call axissymstok_slp_coef_nmode_r64(nc, zc(1:nc), q, Yq, M, mu, &
             C1(1:3*nc,1:3*q,:), C2(1:3*nc,1:3*q,:), C3(1:3*nc,1:3*q,:))
        call sdspecialquad_r64(nc, zc(1:nc), q, Yq, nvq, wxpq, za, zb, iside, &
             As(1:q,1:nc), Ad(1:q,1:nc), A1(1:q,1:nc), A2(1:q,1:nc), A3(1:q,1:nc), A4(1:q,1:nc))
        do md = 0, M
          do e = 1, 9
            ar = (e-1)/3 + 1; b = mod(e-1,3) + 1
            do iq = 1, q
              do ia = 1, nc
                cc1 = C1(3*(ia-1)+ar, 3*(iq-1)+b, md+1)
                cc2 = C2(3*(ia-1)+ar, 3*(iq-1)+b, md+1)
                cc3 = C3(3*(ia-1)+ar, 3*(iq-1)+b, md+1)
                Gcq(ia,iq) = twopi*cc1*As(iq,ia) + cc2*wsq(iq) &
                           + cmplx(twopi*real(Ad(iq,ia)*(cc3/nvq(iq)), r64), 0.0_r64, r64)
              end do
            end do
            Gpc(1:nc,1:p) = matmul(Gcq(1:nc,1:q), IP2)        ! near: 1 interp = Ib upsample fold 2p->p
            do ia = 1, nc
              do l = 1, p                                     ! then the refinement fold (identity for middle panels)
                A((ar-1)*nt+ci(ia), (b-1)*nso+cjo+l, md+1) = &
                  A((ar-1)*nt+ci(ia), (b-1)*nso+cjo+l, md+1) + sum(Gpc(ia,1:p)*IPk(:,l))
              end do
            end do
          end do
        end do
      end if
      ! ---- far: explicit analytic mode-n two-carrier on the refined p nodes, fold via IPk ----
      do i = 1, nt
        if (cl(i)) cycle
        rt = real(tx(i),r64); zti = aimag(tx(i))
        do md = 0, M
          rn = real(md,r64); n2 = rn*rn; fn2 = 4.0_r64*n2 - 1.0_r64
          do jp = 1, p
            rhop = real(Yp(jp),r64); zh = zti - aimag(Yp(jp)); rho = rt
            rr2 = (rho-rhop)**2 + zh*zh; chi = 1.0_r64 + rr2/(2.0_r64*rho*rhop)
            rr = rho*rhop; rrt = (rho*rhop)**(-1.5_r64); srt = 1.0_r64/sqrt(rho*rhop)
            cm1 = chi - 1.0_r64; cp1 = chi + 1.0_r64
            call carrier_r64(chi, md, vk, ve, Fn, An, dFn)
            sk1 = cmplx(-(1.0_r64/8.0_r64)*rhop*ipi*rrt*(rho*rho+rhop*rhop-4.0_r64*chi*rr), 0.0_r64, r64)
            se1 = cmplx(-(1.0_r64/8.0_r64)*rhop*ipi*rrt*(1.0_r64/cm1)*(1.0_r64/fn2)* &
                  ( 2.0_r64*rr + 4.0_r64*n2*rr + chi*rho*rho + chi*rhop*rhop - 4.0_r64*chi*chi*rr &
                  + 4.0_r64*n2*chi*chi*rr - 4.0_r64*chi*n2*rho*rho - 4.0_r64*chi*n2*rhop*rhop ), 0.0_r64, r64)
            sk2 = -(1.0_r64/4.0_r64)*rn*srt*ipi*ic*(rho - chi*rhop)
            se2 = -(3.0_r64/4.0_r64)*rn*rhop*srt*ipi*ic*cp1*(1.0_r64/fn2)
            sk3 = cmplx((1.0_r64/8.0_r64)*rhop*rhop*ipi*rrt*zh, 0.0_r64, r64)
            se3 = cmplx((1.0_r64/8.0_r64)*rhop*ipi*rrt*(1.0_r64/cm1)*zh*(rho - chi*rhop), 0.0_r64, r64)
            sk4 = (1.0_r64/4.0_r64)*rn*rhop*rhop*ipi*rrt*ic*(rhop - chi*rho)
            se4 = (3.0_r64/4.0_r64)*rn*rhop*srt*ipi*ic*cp1*(1.0_r64/fn2)
            sk5 = cmplx((1.0_r64/2.0_r64)*chi*rhop*srt*ipi, 0.0_r64, r64)
            se5 = cmplx(-(1.0_r64/2.0_r64)*rhop*srt*ipi*cp1*(n2-1.0_r64)*(1.0_r64/fn2), 0.0_r64, r64)
            sk6 = -(1.0_r64/4.0_r64)*ic*rn*rhop*rhop*ipi*rrt*zh
            sk7 = cmplx(-(1.0_r64/8.0_r64)*srt*ipi*zh, 0.0_r64, r64)
            se7 = cmplx(-(1.0_r64/8.0_r64)*rhop*ipi*rrt*(1.0_r64/cm1)*zh*(rhop - chi*rho), 0.0_r64, r64)
            sk8 = -(1.0_r64/4.0_r64)*ic*rn*srt*ipi*zh
            sk9 = cmplx((1.0_r64/4.0_r64)*rhop*srt*ipi, 0.0_r64, r64)
            se9 = cmplx((1.0_r64/8.0_r64)*rhop*ipi*rrt*(1.0_r64/cm1)*zh*zh, 0.0_r64, r64)
            ws = wsp(jp)*muinv
            ff(jp,1) = (sk1*vk + se1*ve)*ws; ff(jp,2) = (sk2*vk + se2*ve)*ws; ff(jp,3) = (sk3*vk + se3*ve)*ws
            ff(jp,4) = (sk4*vk + se4*ve)*ws; ff(jp,5) = (sk5*vk + se5*ve)*ws; ff(jp,6) = (sk6*vk)*ws
            ff(jp,7) = (sk7*vk + se7*ve)*ws; ff(jp,8) = (sk8*vk)*ws;          ff(jp,9) = (sk9*vk + se9*ve)*ws
          end do
          do e = 1, 9
            ar = (e-1)/3 + 1; b = mod(e-1,3) + 1
            do l = 1, p
              A((ar-1)*nt+i, (b-1)*nso+cjo+l, md+1) = &
                A((ar-1)*nt+i, (b-1)*nso+cjo+l, md+1) + sum(ff(1:p,e)*IPk(:,l))
            end do
          end do
        end do
      end do
    end do
    deallocate(tin, cl, ci, zc, C1, C2, C3, As, Ad, A1, A2, A3, A4, Gcq, Gpc, ff)
  end subroutine axissymstok_slp_blockmat_nmode_r64

  subroutine axissymstok_slpn_blockmat_r64(nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, mu, A)
    integer(8),   intent(in)    :: nt, p, np, iside, iclosed
    complex(r64), intent(in)    :: tx(nt), tnx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: sws(p*np), tpan(np+1), mu
    real(r64),    intent(inout) :: A(2*nt, 2*np*p)
    integer(8) :: q, nso, i, j, jp, iq, ia, l, nc, cjo, k, jo, ne, npa, nq_az
    real(r64)  :: twopi, sumws, rho, rhop, zh, drho, rr2, chi, ws, spd, rt, zti, nr, nz, vk, ve, Fn, An, dFn
    real(r64)  :: t1, tN, tm, denom, rlo, rhi, tgi, split, thr, vrr, vrz, vzr, vzz
    complex(r64) :: ic, za, zb, Dz
    real(r64)    :: tglp(p), wglp(p), Dp(p,p), tglq(2*p), wglq(2*p), Dq(2*p,2*p), IP2(2*p,p)
    real(r64)    :: tc(p,np), rk(p), re2(2), IPk(p,p), IPe2(2,p), wsq(2*p), wsp(p)
    complex(r64) :: Yp(p), Yq(2*p), dYp(p), dYq(2*p), nvq(2*p), wxpq(2*p), ec2(2), xc(p,np)
    real(r64),    allocatable :: tin(:), tt_az(:), ww_az(:)
    logical,      allocatable :: cl(:)
    integer(8),   allocatable :: ci(:)
    complex(r64), allocatable :: zc(:), zcn(:)
    real(r64),    allocatable :: C1(:,:), C2(:,:), As(:,:), A1(:,:), A2(:,:), A3(:,:), A4(:,:), Gcl(:,:), Gp(:,:), ff(:,:)
    complex(r64), allocatable :: C3(:,:), C4(:,:), Ad(:,:)
    q = 2*p; nso = np*p; ic = (0.0_r64,1.0_r64); twopi = 2.0_r64*acos(-1.0_r64)
    call gauss_r64(p,     tglp, wglp, Dp)
    call gauss_r64(2_8*p, tglq, wglq, Dq)
    call lagrange_interp_r64(p, tglp, q, tglq, IP2)
    do j = 1, np
      do i = 1, p
        xc(i,j) = sx((j-1)*p+i)
        tc(i,j) = tpan(j) + (1.0_r64+tglp(i))/2.0_r64*(tpan(j+1)-tpan(j))
      end do
    end do
    t1 = tc(1,1); tN = tc(p,np)
    ! ---- dyadic pole refinement of tpan -> tin (mirror StoSLPnAxiBlockMat_v0: 1st/last only) ----
    allocate(tin(np+1+256)); ne = np+1; tin(1:ne) = tpan(1:np+1)
    do while (tin(2) > t1)
      split = 0.5_r64*(tin(1)+tin(2))
      do i = ne, 2, -1; tin(i+1) = tin(i); end do
      tin(2) = split; ne = ne+1
    end do
    do while (tin(ne-1) < tN)
      split = 0.5_r64*(tin(ne)+tin(ne-1))
      tin(ne+1) = tin(ne); tin(ne) = split; ne = ne+1
    end do
    npa = ne-1
    thr = 0.1_r64
    call weight_setup_r64(p, 32_8*p, 1.0e-7_r64, nq_az, tt_az, ww_az)   ! graded azimuthal mesh (far close pairs)
    allocate(cl(nt), ci(nt), zc(nt), zcn(nt), C1(2*nt,2*q), C2(2*nt,2*q), C3(2*nt,2*q), C4(2*nt,2*q))
    allocate(As(q,nt), Ad(q,nt), A1(q,nt), A2(q,nt), A3(q,nt), A4(q,nt), Gcl(2*nt,2*q), Gp(nt,p), ff(p,4))
    A = 0.0_r64
    do k = 1, npa
      ! ---- this refined panel: parent jo, coarse->refined interp IPk, positions/endpoints ----
      tm = 0.5_r64*(tin(k)+tin(k+1)); jo = 1
      do j = 1, np
        if (tm >= tpan(j) .and. tm < tpan(j+1)) jo = j
      end do
      denom = tpan(jo+1) - tpan(jo)
      do i = 1, p
        tgi = tin(k) + (1.0_r64+tglp(i))/2.0_r64*(tin(k+1)-tin(k))
        rk(i) = (2.0_r64*tgi - (tpan(jo)+tpan(jo+1)))/denom
      end do
      call lagrange_interp_r64(p, tglp, p, rk, IPk)
      Yp = matmul(IPk, xc(:,jo))
      rlo = (2.0_r64*tin(k)   - (tpan(jo)+tpan(jo+1)))/denom
      rhi = (2.0_r64*tin(k+1) - (tpan(jo)+tpan(jo+1)))/denom
      re2(1) = rlo; re2(2) = rhi
      call lagrange_interp_r64(p, tglp, 2_8, re2, IPe2)
      ec2 = matmul(IPe2, xc(:,jo)); za = ec2(1); zb = ec2(2)
      cjo = (jo-1)*p
      ! ---- refined-panel geometry (upsample + speed/source normal via D-matrix) ----
      Yq = matmul(IP2, Yp); dYq = matmul(Dq, Yq)
      do iq = 1, q
        spd = abs(dYq(iq)); nvq(iq) = -ic*dYq(iq)/spd; wsq(iq) = wglq(iq)*spd; wxpq(iq) = dYq(iq)*wglq(iq)
      end do
      dYp = matmul(Dp, Yp)
      do jp = 1, p
        wsp(jp) = wglp(jp)*abs(dYp(jp))
      end do
      sumws = sum(wsp)
      nc = 0
      do i = 1, nt
        cl(i) = (abs(tx(i)-za) + abs(tx(i)-zb)) < 1.5_r64*sumws
        if (cl(i)) then; nc = nc + 1; ci(nc) = i; zc(nc) = tx(i); zcn(nc) = tnx(i); end if
      end do
      ! ---- near: 2p S' close (coef + sdspecialquad, tgt normal), fold q->p (IP2) -> coarse (IPk) ----
      if (nc > 0) then
        call axissymstok_slpn_coef_r64(nc, zc(1:nc), zcn(1:nc), q, Yq, mu, &
             C1(1:2*nc,1:2*q), C2(1:2*nc,1:2*q), C3(1:2*nc,1:2*q), C4(1:2*nc,1:2*q))
        call sdspecialquad_r64(nc, zc(1:nc), q, Yq, nvq, wxpq, za, zb, iside, &
             As(1:q,1:nc), Ad(1:q,1:nc), A1(1:q,1:nc), A2(1:q,1:nc), A3(1:q,1:nc), A4(1:q,1:nc))
        do iq = 1, q
          do ia = 1, nc
            Dz = cmplx(A1(iq,ia), -A2(iq,ia), r64)
            Gcl(ia,    iq)   = twopi*C1(ia,iq)*As(iq,ia)     + C2(ia,iq)*wsq(iq) &
                 + twopi*real(Ad(iq,ia)*(C3(ia,iq)/nvq(iq)), r64) - twopi*real(Dz*(C4(ia,iq)/nvq(iq)), r64)
            Gcl(ia,    q+iq) = twopi*C1(ia,q+iq)*As(iq,ia)   + C2(ia,q+iq)*wsq(iq) &
                 + twopi*real(Ad(iq,ia)*(C3(ia,q+iq)/nvq(iq)), r64) - twopi*real(Dz*(C4(ia,q+iq)/nvq(iq)), r64)
            Gcl(nc+ia, iq)   = twopi*C1(nc+ia,iq)*As(iq,ia)  + C2(nc+ia,iq)*wsq(iq) &
                 + twopi*real(Ad(iq,ia)*(C3(nc+ia,iq)/nvq(iq)), r64) - twopi*real(Dz*(C4(nc+ia,iq)/nvq(iq)), r64)
            Gcl(nc+ia, q+iq) = twopi*C1(nc+ia,q+iq)*As(iq,ia)+ C2(nc+ia,q+iq)*wsq(iq) &
                 + twopi*real(Ad(iq,ia)*(C3(nc+ia,q+iq)/nvq(iq)), r64) - twopi*real(Dz*(C4(nc+ia,q+iq)/nvq(iq)), r64)
          end do
        end do
        Gp(1:nc,:) = matmul(Gcl(1:nc, 1:q), IP2)
        do ia = 1, nc; do l = 1, p; A(ci(ia),    cjo+l)     = A(ci(ia),    cjo+l)     + sum(Gp(ia,1:p)*IPk(:,l)); end do; end do
        Gp(1:nc,:) = matmul(Gcl(1:nc, q+1:2*q), IP2)
        do ia = 1, nc; do l = 1, p; A(ci(ia),    nso+cjo+l) = A(ci(ia),    nso+cjo+l) + sum(Gp(ia,1:p)*IPk(:,l)); end do; end do
        Gp(1:nc,:) = matmul(Gcl(nc+1:2*nc, 1:q), IP2)
        do ia = 1, nc; do l = 1, p; A(nt+ci(ia), cjo+l)     = A(nt+ci(ia), cjo+l)     + sum(Gp(ia,1:p)*IPk(:,l)); end do; end do
        Gp(1:nc,:) = matmul(Gcl(nc+1:2*nc, q+1:2*q), IP2)
        do ia = 1, nc; do l = 1, p; A(nt+ci(ia), nso+cjo+l) = A(nt+ci(ia), nso+cjo+l) + sum(Gp(ia,1:p)*IPk(:,l)); end do; end do
      end if
      ! ---- far: slpn9 meridian (negated) on the refined p nodes, tgt normal, fold via IPk ----
      do i = 1, nt
        if (cl(i)) cycle
        rt = real(tx(i),r64); zti = aimag(tx(i)); nr = real(tnx(i),r64); nz = aimag(tnx(i))
        do jp = 1, p
          rhop = real(Yp(jp),r64); zh = zti - aimag(Yp(jp))
          rho = rt; drho = rho - rhop; rr2 = drho*drho + zh*zh
          if (sqrt(rr2) < thr) then
            ! close pair in the far bucket: cancellation-free azimuthal traction quad
            call axissymstok_slpn_aziquad_r64(rho, rhop, zh, nr, nz, nq_az, tt_az, ww_az, vrr, vrz, vzr, vzz)
            ff(jp,1) = vrr*rhop*wsp(jp); ff(jp,2) = vrz*rhop*wsp(jp)
            ff(jp,3) = vzr*rhop*wsp(jp); ff(jp,4) = vzz*rhop*wsp(jp)
            cycle
          end if
          chi = 1.0_r64 + rr2/(2.0_r64*rho*rhop)
          call carrier_r64(chi, 0_8, vk, ve, Fn, An, dFn)
          block
            real(r64) :: t2,t3,t4,t5,t6,t7,t8,t9,t11,t13,t15,t16,t17,t18,t19,t20,t23,t24,t25,t26,t29,t30,t31
            real(r64) :: t32,t33,t38,t39,t44,t45,t46,t47,t52,t55,t56,t59,t62,t66
            real(r64) :: c_1,c_2,c_11,c_12,c_31,c_32,c_41,c_42
            t2=nr*rho; t3=nr*rhop; t4=rho*rhop; t5=nz*zh; t6=chi+1.0_r64; t7=chi*chi; t8=chi**3
            t11=rho*rho; t13=rhop*rhop; t15=zh*zh; t17=1.0_r64/acos(-1.0_r64); t20=chi-1.0_r64
            t9=t7*t7; t16=t3*3.0_r64; t18=rho*t5; t19=rhop*t5; t23=-t3; t24=rho*t2; t25=t2*t11
            t26=rhop*t3; t29=chi*rhop*t2*2.0_r64; t30=t5*t11; t31=1.0_r64/t6
            t32=t7-1.0_r64; t33=t4*t5*4.0_r64; t44=t2*t13*3.0_r64; t45=t2*t4*7.0_r64; t46=1.0_r64/t20; t52=t2+t5
            t55=rhop*t2*t7*2.0_r64; t59=1.0_r64/t4**2.5_r64; t38=t24*3.0_r64; t39=rhop*t16; t47=t46*t46
            t56=-1.0_r64; t62=1.0_r64/t32; t66=1.0_r64/t56
            c_1  = rhop*t17*t59*t62*(-t33-t45+chi*t25+chi*t30-t3*t13*3.0_r64+chi*t2*t13*11.0_r64+chi*t5*t13 &
                 +t2*t4*t7*4.0_r64+t4*t5*t7*2.0_r64-t2*t8*t13*8.0_r64+t3*t7*t13*2.0_r64)*(-1.0_r64/8.0_r64)
            c_2  = rhop*t17*t31*t47*t59*t66*(t25*3.0_r64+t30*3.0_r64+t44+t5*t13*3.0_r64+t7*t25 &
                 +t7*t30-chi*t2*t4*16.0_r64-chi*t4*t5*10.0_r64-chi*t3*t13*6.0_r64+t2*t4*t8*4.0_r64 &
                 +t4*t5*t8*2.0_r64+t2*t7*t13*17.0_r64-t2*t9*t13*8.0_r64+t3*t8*t13*2.0_r64 &
                 +t5*t7*t13)*(-1.0_r64/8.0_r64)
            c_11 = (rhop*t17*t59*t62*zh*(-t18-t24-t26*3.0_r64+t29+chi*t19+t7*t26*2.0_r64))/8.0_r64
            c_12 = rhop*t17*t31*t47*t59*zh*(t19*3.0_r64+t55-chi*t18*4.0_r64-chi*t24*4.0_r64-chi*t26*6.0_r64 &
                 +rhop*t2*6.0_r64+t7*t19+t8*t26*2.0_r64)*(-1.0_r64/8.0_r64)
            c_31 = rhop*t17*t59*t62*zh*(-t19+t55+chi*t18+chi*t24+chi*t26-rhop*t2*4.0_r64)*(-1.0_r64/8.0_r64)
            c_32 = (rhop*t17*t31*t47*t59*zh*(t18*3.0_r64+t38+t39-chi*t19*4.0_r64+t7*t18+t7*t24+t7*t26 &
                 -chi*rhop*t2*10.0_r64+rhop*t2*t8*2.0_r64))/8.0_r64
            c_41 = (rhop*t15*t17*t56*t59*t62*(t52+chi*t23))/8.0_r64
            c_42 = rhop*t15*t17*t31*t47*t59*(t16-chi*t2*4.0_r64-chi*t5*4.0_r64+t3*t7)*(-1.0_r64/8.0_r64)
            ff(jp,1) = -(c_1 *vk + c_2 *ve)*wsp(jp); ff(jp,2) = -(c_11*vk + c_12*ve)*wsp(jp)
            ff(jp,3) = -(c_31*vk + c_32*ve)*wsp(jp); ff(jp,4) = -(c_41*vk + c_42*ve)*wsp(jp)
          end block
        end do
        do l = 1, p
          A(i,     cjo+l)     = A(i,     cjo+l)     + sum(ff(1:p,1)*IPk(:,l))
          A(i,     nso+cjo+l) = A(i,     nso+cjo+l) + sum(ff(1:p,2)*IPk(:,l))
          A(nt+i,  cjo+l)     = A(nt+i,  cjo+l)     + sum(ff(1:p,3)*IPk(:,l))
          A(nt+i,  nso+cjo+l) = A(nt+i,  nso+cjo+l) + sum(ff(1:p,4)*IPk(:,l))
        end do
      end do
    end do
    deallocate(tin, tt_az, ww_az, cl, ci, zc, zcn, C1, C2, C3, C4, As, Ad, A1, A2, A3, A4, Gcl, Gp, ff)
  end subroutine axissymstok_slpn_blockmat_r64

  subroutine axissymstok_slpn_blockmat_nmode_r64(nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, mu, A)
    integer(8),   intent(in)    :: nt, p, np, M, iside, iclosed
    complex(r64), intent(in)    :: tx(nt), tnx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: sws(p*np), tpan(np+1), mu
    complex(r64), intent(inout) :: A(3*nt, 3*np*p, M+1)
    integer(8), parameter :: Ksub = 4
    integer(8) :: q, pq, nso, i, j, jp, iq, ia, l, e, md, ar, b, c, qo, ne, npa, nk, nkc, jj, cols
    real(r64)  :: twopi, sumwso, sumwsc, rho, rhop, zh, rr2, chi, ws, spd, rt, zti, nr, nz, vk, ve, Fn, An, dFn, rn
    real(r64)  :: st1, stN, tm, denom, rlo, rhi, tgi, split
    complex(r64) :: ic, zac, zbc, cc1, cc2, cc3, cc4, Slog, Dval, Dz, Acau, Ahy
    real(r64)    :: tglp(p), wglp(p), Dp(p,p), tglq(2*p), wglq(2*p), Dq(2*p,2*p), IP2(2*p,p)
    real(r64)    :: tc(p,np), rk(p), re2(2), Lc(p,p), IPe2(2,p), IPqc(2*p,p), wsq(2*p), wsp(p)
    complex(r64) :: Ypb(p), Yq(2*p), dYq(2*p), dYp(p), nvq(2*p), wxpq(2*p), ec2(2), xc(p,np)
    real(r64),    allocatable :: tin(:)
    integer(8),   allocatable :: joa(:), ci(:), nidx(:)
    logical,      allocatable :: ikq(:), ikc(:)
    complex(r64), allocatable :: zc(:), zcn(:), znear(:), znearn(:), C1(:,:,:), C2(:,:,:), C3(:,:,:), C4(:,:,:), Gcq(:,:), Gpc(:,:), ff(:,:)
    real(r64),    allocatable :: As(:,:), A1(:,:), A2(:,:), A3(:,:), A4(:,:)
    complex(r64), allocatable :: Ad(:,:)
    q = 2*p; nso = np*p; ic = (0.0_r64,1.0_r64); twopi = 2.0_r64*acos(-1.0_r64)
    call gauss_r64(p,     tglp, wglp, Dp)
    call gauss_r64(2_8*p, tglq, wglq, Dq)
    call lagrange_interp_r64(p, tglp, q, tglq, IP2)
    do j = 1, np
      do i = 1, p
        xc(i,j) = sx((j-1)*p+i)
        tc(i,j) = tpan(j) + (1.0_r64+tglp(i))/2.0_r64*(tpan(j+1)-tpan(j))
      end do
    end do
    st1 = tc(1,1); stN = tc(p,np)
    ! ---- two-level aux mesh: Ksub uniform subpanels per panel + dyadic pole refine ----
    allocate(tin(np*Ksub+1+512)); ne = 0
    do j = 1, np
      do i = 0, Ksub-1
        ne = ne+1; tin(ne) = tpan(j) + (tpan(j+1)-tpan(j))*real(i,r64)/real(Ksub,r64)
      end do
    end do
    ne = ne+1; tin(ne) = tpan(np+1)
    do while (tin(2) > st1)
      split = 0.5_r64*(tin(1)+tin(2)); do i = ne, 2, -1; tin(i+1) = tin(i); end do; tin(2) = split; ne = ne+1
    end do
    do while (tin(ne-1) < stN)
      split = 0.5_r64*(tin(ne)+tin(ne-1)); tin(ne+1) = tin(ne); tin(ne) = split; ne = ne+1
    end do
    npa = ne-1
    allocate(joa(npa))
    do c = 1, npa
      tm = 0.5_r64*(tin(c)+tin(c+1)); qo = 1
      do j = 1, np
        if (tm >= tpan(j) .and. tm < tpan(j+1)) qo = j
      end do
      joa(c) = qo
    end do
    allocate(ci(nt), zc(nt), zcn(nt), nidx(nt), znear(nt), znearn(nt), ikq(nt), ikc(nt))
    allocate(C1(3*nt,3*q,M+1), C2(3*nt,3*q,M+1), C3(3*nt,3*q,M+1), C4(3*nt,3*q,M+1))
    allocate(As(q,nt), Ad(q,nt), A1(q,nt), A2(q,nt), A3(q,nt), A4(q,nt), Gcq(nt,q), Gpc(nt,p), ff(p,9))
    A = (0.0_r64,0.0_r64)
    do pq = 1, np                                                   ! ===== outer level: original panel pq =====
      cols = (pq-1)*p
      sumwso = 0.0_r64
      do jp = 1, p
        sumwso = sumwso + sws(cols+jp)
      end do
      nk = 0
      do i = 1, nt
        ikq(i) = (abs(tx(i)-sxlo(pq)) + abs(tx(i)-sxhi(pq))) < 1.85_r64*sumwso
        if (ikq(i)) then; nk = nk + 1; nidx(nk) = i; znear(nk) = tx(i); znearn(nk) = tnx(i); end if
      end do
      ! ---- outer far: naive slpn9far (negated) on the COARSE original nodes, TARGET normals ----
      do i = 1, nt
        if (ikq(i)) cycle
        rt = real(tx(i),r64); zti = aimag(tx(i)); nr = real(tnx(i),r64); nz = aimag(tnx(i))
        do jp = 1, p
          jj = cols + jp; rhop = real(sx(jj),r64); zh = zti - aimag(sx(jj)); ws = sws(jj)
          rho = rt; rr2 = (rho-rhop)**2 + zh*zh; chi = 1.0_r64 + rr2/(2.0_r64*rho*rhop)
          do md = 0, M
            rn = real(md,r64)
            call carrier_r64(chi, md, vk, ve, Fn, An, dFn)
            block
              real(r64) :: t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t13,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25
              real(r64) :: t26,t27,t28,t29,t30,t31,t32,t33,t38,t39,t41,t42,t43,t44,t45,t46,t47,t52,t53,t55,t56,t59,t61,t62,t66
              complex(r64) :: t34,t35,t36,t37,t48,t50,t51,t54,t57,t58,t60,t63,t64,t65,t69, coK(9), coE(9)
              t2=nr*rho; t3=nr*rhop; t4=rho*rhop; t5=nz*zh; t6=chi+1.0_r64; t7=chi*chi; t8=chi**3
              t10=rn*rn; t11=rho*rho; t13=rhop*rhop; t15=zh*zh; t17=1.0_r64/acos(-1.0_r64); t20=chi-1.0_r64; t21=-rhop
              t9=t7*t7; t16=t3*3.0_r64; t18=rho*t5; t19=rhop*t5; t22=t10*4.0_r64; t23=-t3; t24=rho*t2; t25=t2*t11
              t26=rhop*t3; t27=chi*t3*4.0_r64; t28=rhop*t2*4.0_r64; t29=chi*rhop*t2*2.0_r64; t30=t5*t11; t31=1.0_r64/t6
              t32=t7-1.0_r64; t33=t4*t5*4.0_r64; t34=t2*ic; t35=t2*2.0_r64*ic; t36=t3*ic; t37=t5*ic; t41=rho+t21
              t44=t2*t13*3.0_r64; t45=t2*t4*7.0_r64; t46=1.0_r64/t20; t50=chi*t3*2.0_r64*ic; t52=t2+t5
              t55=rhop*t2*t7*2.0_r64; t59=1.0_r64/t4**2.5_r64; t60=chi*t3*(-ic); t63=rhop*t2*t10*4.0_r64*ic
              t38=t24*3.0_r64; t39=rhop*t16; t42=-t27; t43=-t28; t47=t46*t46; t48=-t36; t51=t18*ic; t53=t41*t41
              t54=chi*rhop*t35; t56=t22-1.0_r64; t57=t24*ic; t58=t26*ic; t61=1.0_r64/t41; t62=1.0_r64/t32
              t64=t10*t18*4.0_r64*ic; t65=t10*t24*4.0_r64*ic; t69=t34+t37+t60; t66=1.0_r64/t56
              coK(1) = cmplx(rhop*t17*t59*t62*(-t33-t45+chi*t25+chi*t30-t3*t13*3.0_r64+chi*t2*t13*11.0_r64+chi*t5*t13 &
                     -chi*t10*t25*4.0_r64-chi*t10*t30*4.0_r64+t2*t4*t7*4.0_r64+t4*t5*t7*2.0_r64-t2*t8*t13*8.0_r64 &
                     +t3*t7*t13*2.0_r64+t2*t4*t22+t4*t5*t22-chi*t2*t10*t13*8.0_r64-chi*t5*t10*t13*4.0_r64 &
                     +t2*t4*t7*t10*8.0_r64-t2*t8*t10*t13*4.0_r64+t4*t5*t7*t22+t3*t7*t13*t22)*(-1.0_r64/8.0_r64), 0.0_r64, r64)
              coE(1) = cmplx(rhop*t17*t31*t47*t59*t66*(t25*3.0_r64+t30*3.0_r64+t44+t5*t13*3.0_r64+t7*t25-t10*t25*12.0_r64 &
                     +t7*t30-t10*t30*12.0_r64-chi*t2*t4*16.0_r64-chi*t4*t5*10.0_r64-chi*t3*t13*6.0_r64+t2*t4*t8*4.0_r64 &
                     +t4*t5*t8*2.0_r64+t2*t7*t13*17.0_r64-t2*t9*t13*8.0_r64+t3*t8*t13*2.0_r64-t2*t10*t13*24.0_r64 &
                     +t5*t7*t13-t5*t10*t13*12.0_r64-t7*t10*t25*4.0_r64-t7*t10*t30*4.0_r64+chi*t2*t4*t10*64.0_r64 &
                     +chi*t4*t5*t10*40.0_r64+chi*t3*t10*t13*24.0_r64-t2*t4*t8*t10*16.0_r64-t4*t5*t8*t10*8.0_r64 &
                     -t2*t7*t10*t13*44.0_r64+t2*t9*t10*t13*20.0_r64-t3*t8*t10*t13*8.0_r64-t5*t7*t10*t13*4.0_r64)*(-1.0_r64/8.0_r64), 0.0_r64, r64)
              coK(2) = rn*rhop/t4**1.5_r64*t17*(t35+t37-t50)*(-3.0_r64/4.0_r64)
              coE(2) = rn/t4**1.5_r64*t17*t46*t66*(t26*(-3.0_r64*ic)-t51+t54-t57+t64+t65+chi*t19*ic+t7*t26*2.0_r64*ic &
                     -chi*t10*t19*4.0_r64*ic+t7*t10*t26*4.0_r64*ic-chi*rhop*t2*t10*8.0_r64*ic)*(-1.0_r64/4.0_r64)
              coK(3) = cmplx((rhop*t17*t59*t62*zh*(-t18-t24-t26*3.0_r64+t29+chi*t19+t7*t26*2.0_r64+t18*t22+t22*t24 &
                     -chi*t10*t19*4.0_r64+t7*t22*t26-chi*rhop*t2*t10*8.0_r64))/8.0_r64, 0.0_r64, r64)
              coE(3) = cmplx(rhop*t17*t31*t47*t59*zh*(t19*3.0_r64+t55-chi*t18*4.0_r64-chi*t24*4.0_r64-chi*t26*6.0_r64 &
                     +rhop*t2*6.0_r64+t7*t19+t8*t26*2.0_r64)*(-1.0_r64/8.0_r64), 0.0_r64, r64)
              coK(4) = rn*t13*t17*t59*(t51+t57+t58-chi*rhop*t2*2.0_r64*ic)*(3.0_r64/4.0_r64)
              coE(4) = (rn*t13*t17*t46*t59*t66*(t19*(-ic)+t63+chi*t51+chi*t57+chi*t58-rhop*t2*4.0_r64*ic+t10*t19*4.0_r64*ic &
                     +t7*t63-chi*t10*t18*4.0_r64*ic-chi*t10*t24*4.0_r64*ic-chi*t10*t26*4.0_r64*ic+rhop*t7*t35))/4.0_r64
              coK(5) = cmplx(rhop/t4**1.5_r64*t17*(t42+t52+t2*t10*2.0_r64+t5*t10*2.0_r64-chi*t3*t10*2.0_r64)*(-1.0_r64/4.0_r64), 0.0_r64, r64)
              coE(5) = cmplx(rhop/t4**1.5_r64*t17*t46*t66*(t16+chi*t2+chi*t5-t3*t7*4.0_r64-t3*t10*6.0_r64 &
                     -chi*t2*t10*4.0_r64-chi*t5*t10*4.0_r64+t3*t7*t10*10.0_r64)*(-1.0_r64/4.0_r64), 0.0_r64, r64)
              coK(6) = rn*t3*t13*t17*t59*zh*cmplx(0.0_r64,-0.75_r64,r64)
              coE(6) = rn*t13*t17*t46*t59*t69*zh*(-1.0_r64/4.0_r64)
              coK(7) = cmplx(rhop*t17*t59*t62*zh*(-t19+t43+t55+chi*t18+chi*t24+chi*t26+t19*t22-chi*t10*t18*4.0_r64 &
                     -chi*t10*t24*4.0_r64-chi*t10*t26*4.0_r64+rhop*t2*t22+rhop*t2*t7*t22)*(-1.0_r64/8.0_r64), 0.0_r64, r64)
              coE(7) = cmplx((rhop*t17*t31*t47*t59*zh*(t18*3.0_r64+t38+t39-chi*t19*4.0_r64+t7*t18+t7*t24+t7*t26 &
                     -chi*rhop*t2*10.0_r64+rhop*t2*t8*2.0_r64))/8.0_r64, 0.0_r64, r64)
              coK(8) = rn*t2*t13*t17*t59*zh*cmplx(0.0_r64,-0.75_r64,r64)
              coE(8) = rn/t4**1.5_r64*t17*t46*t69*zh*(-1.0_r64/4.0_r64)
              coK(9) = cmplx((rhop*t15*t17*t56*t59*t62*(t52+chi*t23))/8.0_r64, 0.0_r64, r64)
              coE(9) = cmplx(rhop*t15*t17*t31*t47*t59*(t16-chi*t2*4.0_r64-chi*t5*4.0_r64+t3*t7)*(-1.0_r64/8.0_r64), 0.0_r64, r64)
              do e = 1, 9
                ar = (e-1)/3 + 1; b = mod(e-1,3) + 1
                A((ar-1)*nt+i, (b-1)*nso+jj, md+1) = -(coK(e)*vk + coE(e)*ve)*ws
              end do
            end block
          end do
        end do
      end do
      if (nk == 0) cycle
      ! ---- outer near: recurse into the Ksub(+dyadic) subpanels of pq, inner re-detect ----
      do c = 1, npa
        if (joa(c) /= pq) cycle
        denom = tpan(pq+1) - tpan(pq)
        do i = 1, p
          tgi = tin(c) + (1.0_r64+tglp(i))/2.0_r64*(tin(c+1)-tin(c))
          rk(i) = (2.0_r64*tgi - (tpan(pq)+tpan(pq+1)))/denom
        end do
        call lagrange_interp_r64(p, tglp, p, rk, Lc)
        Ypb = matmul(Lc, xc(:,pq))
        rlo = (2.0_r64*tin(c)   - (tpan(pq)+tpan(pq+1)))/denom
        rhi = (2.0_r64*tin(c+1) - (tpan(pq)+tpan(pq+1)))/denom
        re2(1) = rlo; re2(2) = rhi
        call lagrange_interp_r64(p, tglp, 2_8, re2, IPe2)
        ec2 = matmul(IPe2, xc(:,pq)); zac = ec2(1); zbc = ec2(2)
        IPqc = matmul(IP2, Lc)
        Yq = matmul(IP2, Ypb); dYq = matmul(Dq, Yq)
        do iq = 1, q
          spd = abs(dYq(iq)); nvq(iq) = -ic*dYq(iq)/spd; wsq(iq) = wglq(iq)*spd; wxpq(iq) = dYq(iq)*wglq(iq)
        end do
        dYp = matmul(Dp, Ypb)
        do jp = 1, p
          wsp(jp) = wglp(jp)*abs(dYp(jp))
        end do
        sumwsc = sum(wsp)
        nkc = 0
        do i = 1, nk
          ikc(i) = (abs(znear(i)-zac) + abs(znear(i)-zbc)) < 1.5_r64*sumwsc
          if (ikc(i)) then; nkc = nkc + 1; ci(nkc) = i; zc(nkc) = znear(i); zcn(nkc) = znearn(i); end if
        end do
        ! ---- inner near: 2p S' close (coef + sdspecialquad), fold q->p_sub (IP2) -> p_coarse (Lc) ----
        if (nkc > 0) then
          call axissymstok_slpn_coef_nmode_r64(nkc, zc(1:nkc), zcn(1:nkc), q, Yq, M, mu, &
               C1(1:3*nkc,1:3*q,:), C2(1:3*nkc,1:3*q,:), C3(1:3*nkc,1:3*q,:), C4(1:3*nkc,1:3*q,:))
          call sdspecialquad_r64(nkc, zc(1:nkc), q, Yq, nvq, wxpq, zac, zbc, iside, &
               As(1:q,1:nkc), Ad(1:q,1:nkc), A1(1:q,1:nkc), A2(1:q,1:nkc), A3(1:q,1:nkc), A4(1:q,1:nkc))
          do md = 0, M
            do e = 1, 9
              ar = (e-1)/3 + 1; b = mod(e-1,3) + 1
              do iq = 1, q
                do ia = 1, nkc
                  cc1 = C1(3*(ia-1)+ar, 3*(iq-1)+b, md+1); cc2 = C2(3*(ia-1)+ar, 3*(iq-1)+b, md+1)
                  cc3 = C3(3*(ia-1)+ar, 3*(iq-1)+b, md+1); cc4 = C4(3*(ia-1)+ar, 3*(iq-1)+b, md+1)
                  Slog = cmplx(As(iq,ia),0.0_r64,r64); Dval = Ad(iq,ia); Dz = cmplx(A1(iq,ia),-A2(iq,ia),r64)
                  if (mod(e,2) == 0) then
                    Acau =  ic*twopi*aimag(Dval*(cc3/nvq(iq)))
                    Ahy  = -ic*twopi*aimag(Dz  *(cc4/nvq(iq)))
                  else
                    Acau = cmplx(twopi*real(Dval*(cc3/nvq(iq)), r64), 0.0_r64, r64)
                    Ahy  = cmplx(-twopi*real(Dz *(cc4/nvq(iq)), r64), 0.0_r64, r64)
                  end if
                  Gcq(ia,iq) = twopi*cc1*Slog + cc2*wsq(iq) + Acau + Ahy
                end do
              end do
              Gpc(1:nkc,1:p) = matmul(Gcq(1:nkc,1:q), IPqc)
              do ia = 1, nkc
                do l = 1, p
                  A((ar-1)*nt+nidx(ci(ia)), (b-1)*nso+cols+l, md+1) = &
                    A((ar-1)*nt+nidx(ci(ia)), (b-1)*nso+cols+l, md+1) + Gpc(ia,l)
                end do
              end do
            end do
          end do
        end if
        ! ---- inner far: naive slpn9far (negated) on the FINER subpanel nodes, TARGET normals, fold via Lc ----
        do i = 1, nk
          if (ikc(i)) cycle
          rt = real(znear(i),r64); zti = aimag(znear(i)); nr = real(znearn(i),r64); nz = aimag(znearn(i))
          do md = 0, M
            rn = real(md,r64)
            do jp = 1, p
              rhop = real(Ypb(jp),r64); zh = zti - aimag(Ypb(jp)); rho = rt
              rr2 = (rho-rhop)**2 + zh*zh; chi = 1.0_r64 + rr2/(2.0_r64*rho*rhop)
              call carrier_r64(chi, md, vk, ve, Fn, An, dFn)
              block
                real(r64) :: t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t13,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25
                real(r64) :: t26,t27,t28,t29,t30,t31,t32,t33,t38,t39,t41,t42,t43,t44,t45,t46,t47,t52,t53,t55,t56,t59,t61,t62,t66
                complex(r64) :: t34,t35,t36,t37,t48,t50,t51,t54,t57,t58,t60,t63,t64,t65,t69, coK(9), coE(9)
                t2=nr*rho; t3=nr*rhop; t4=rho*rhop; t5=nz*zh; t6=chi+1.0_r64; t7=chi*chi; t8=chi**3
                t10=rn*rn; t11=rho*rho; t13=rhop*rhop; t15=zh*zh; t17=1.0_r64/acos(-1.0_r64); t20=chi-1.0_r64; t21=-rhop
                t9=t7*t7; t16=t3*3.0_r64; t18=rho*t5; t19=rhop*t5; t22=t10*4.0_r64; t23=-t3; t24=rho*t2; t25=t2*t11
                t26=rhop*t3; t27=chi*t3*4.0_r64; t28=rhop*t2*4.0_r64; t29=chi*rhop*t2*2.0_r64; t30=t5*t11; t31=1.0_r64/t6
                t32=t7-1.0_r64; t33=t4*t5*4.0_r64; t34=t2*ic; t35=t2*2.0_r64*ic; t36=t3*ic; t37=t5*ic; t41=rho+t21
                t44=t2*t13*3.0_r64; t45=t2*t4*7.0_r64; t46=1.0_r64/t20; t50=chi*t3*2.0_r64*ic; t52=t2+t5
                t55=rhop*t2*t7*2.0_r64; t59=1.0_r64/t4**2.5_r64; t60=chi*t3*(-ic); t63=rhop*t2*t10*4.0_r64*ic
                t38=t24*3.0_r64; t39=rhop*t16; t42=-t27; t43=-t28; t47=t46*t46; t48=-t36; t51=t18*ic; t53=t41*t41
                t54=chi*rhop*t35; t56=t22-1.0_r64; t57=t24*ic; t58=t26*ic; t61=1.0_r64/t41; t62=1.0_r64/t32
                t64=t10*t18*4.0_r64*ic; t65=t10*t24*4.0_r64*ic; t69=t34+t37+t60; t66=1.0_r64/t56
                coK(1) = cmplx(rhop*t17*t59*t62*(-t33-t45+chi*t25+chi*t30-t3*t13*3.0_r64+chi*t2*t13*11.0_r64+chi*t5*t13 &
                       -chi*t10*t25*4.0_r64-chi*t10*t30*4.0_r64+t2*t4*t7*4.0_r64+t4*t5*t7*2.0_r64-t2*t8*t13*8.0_r64 &
                       +t3*t7*t13*2.0_r64+t2*t4*t22+t4*t5*t22-chi*t2*t10*t13*8.0_r64-chi*t5*t10*t13*4.0_r64 &
                       +t2*t4*t7*t10*8.0_r64-t2*t8*t10*t13*4.0_r64+t4*t5*t7*t22+t3*t7*t13*t22)*(-1.0_r64/8.0_r64), 0.0_r64, r64)
                coE(1) = cmplx(rhop*t17*t31*t47*t59*t66*(t25*3.0_r64+t30*3.0_r64+t44+t5*t13*3.0_r64+t7*t25-t10*t25*12.0_r64 &
                       +t7*t30-t10*t30*12.0_r64-chi*t2*t4*16.0_r64-chi*t4*t5*10.0_r64-chi*t3*t13*6.0_r64+t2*t4*t8*4.0_r64 &
                       +t4*t5*t8*2.0_r64+t2*t7*t13*17.0_r64-t2*t9*t13*8.0_r64+t3*t8*t13*2.0_r64-t2*t10*t13*24.0_r64 &
                       +t5*t7*t13-t5*t10*t13*12.0_r64-t7*t10*t25*4.0_r64-t7*t10*t30*4.0_r64+chi*t2*t4*t10*64.0_r64 &
                       +chi*t4*t5*t10*40.0_r64+chi*t3*t10*t13*24.0_r64-t2*t4*t8*t10*16.0_r64-t4*t5*t8*t10*8.0_r64 &
                       -t2*t7*t10*t13*44.0_r64+t2*t9*t10*t13*20.0_r64-t3*t8*t10*t13*8.0_r64-t5*t7*t10*t13*4.0_r64)*(-1.0_r64/8.0_r64), 0.0_r64, r64)
                coK(2) = rn*rhop/t4**1.5_r64*t17*(t35+t37-t50)*(-3.0_r64/4.0_r64)
                coE(2) = rn/t4**1.5_r64*t17*t46*t66*(t26*(-3.0_r64*ic)-t51+t54-t57+t64+t65+chi*t19*ic+t7*t26*2.0_r64*ic &
                       -chi*t10*t19*4.0_r64*ic+t7*t10*t26*4.0_r64*ic-chi*rhop*t2*t10*8.0_r64*ic)*(-1.0_r64/4.0_r64)
                coK(3) = cmplx((rhop*t17*t59*t62*zh*(-t18-t24-t26*3.0_r64+t29+chi*t19+t7*t26*2.0_r64+t18*t22+t22*t24 &
                       -chi*t10*t19*4.0_r64+t7*t22*t26-chi*rhop*t2*t10*8.0_r64))/8.0_r64, 0.0_r64, r64)
                coE(3) = cmplx(rhop*t17*t31*t47*t59*zh*(t19*3.0_r64+t55-chi*t18*4.0_r64-chi*t24*4.0_r64-chi*t26*6.0_r64 &
                       +rhop*t2*6.0_r64+t7*t19+t8*t26*2.0_r64)*(-1.0_r64/8.0_r64), 0.0_r64, r64)
                coK(4) = rn*t13*t17*t59*(t51+t57+t58-chi*rhop*t2*2.0_r64*ic)*(3.0_r64/4.0_r64)
                coE(4) = (rn*t13*t17*t46*t59*t66*(t19*(-ic)+t63+chi*t51+chi*t57+chi*t58-rhop*t2*4.0_r64*ic+t10*t19*4.0_r64*ic &
                       +t7*t63-chi*t10*t18*4.0_r64*ic-chi*t10*t24*4.0_r64*ic-chi*t10*t26*4.0_r64*ic+rhop*t7*t35))/4.0_r64
                coK(5) = cmplx(rhop/t4**1.5_r64*t17*(t42+t52+t2*t10*2.0_r64+t5*t10*2.0_r64-chi*t3*t10*2.0_r64)*(-1.0_r64/4.0_r64), 0.0_r64, r64)
                coE(5) = cmplx(rhop/t4**1.5_r64*t17*t46*t66*(t16+chi*t2+chi*t5-t3*t7*4.0_r64-t3*t10*6.0_r64 &
                       -chi*t2*t10*4.0_r64-chi*t5*t10*4.0_r64+t3*t7*t10*10.0_r64)*(-1.0_r64/4.0_r64), 0.0_r64, r64)
                coK(6) = rn*t3*t13*t17*t59*zh*cmplx(0.0_r64,-0.75_r64,r64)
                coE(6) = rn*t13*t17*t46*t59*t69*zh*(-1.0_r64/4.0_r64)
                coK(7) = cmplx(rhop*t17*t59*t62*zh*(-t19+t43+t55+chi*t18+chi*t24+chi*t26+t19*t22-chi*t10*t18*4.0_r64 &
                       -chi*t10*t24*4.0_r64-chi*t10*t26*4.0_r64+rhop*t2*t22+rhop*t2*t7*t22)*(-1.0_r64/8.0_r64), 0.0_r64, r64)
                coE(7) = cmplx((rhop*t17*t31*t47*t59*zh*(t18*3.0_r64+t38+t39-chi*t19*4.0_r64+t7*t18+t7*t24+t7*t26 &
                       -chi*rhop*t2*10.0_r64+rhop*t2*t8*2.0_r64))/8.0_r64, 0.0_r64, r64)
                coK(8) = rn*t2*t13*t17*t59*zh*cmplx(0.0_r64,-0.75_r64,r64)
                coE(8) = rn/t4**1.5_r64*t17*t46*t69*zh*(-1.0_r64/4.0_r64)
                coK(9) = cmplx((rhop*t15*t17*t56*t59*t62*(t52+chi*t23))/8.0_r64, 0.0_r64, r64)
                coE(9) = cmplx(rhop*t15*t17*t31*t47*t59*(t16-chi*t2*4.0_r64-chi*t5*4.0_r64+t3*t7)*(-1.0_r64/8.0_r64), 0.0_r64, r64)
                do e = 1, 9
                  ff(jp,e) = -(coK(e)*vk + coE(e)*ve)*wsp(jp)
                end do
              end block
            end do
            do e = 1, 9
              ar = (e-1)/3 + 1; b = mod(e-1,3) + 1
              do l = 1, p
                A((ar-1)*nt+nidx(i), (b-1)*nso+cols+l, md+1) = &
                  A((ar-1)*nt+nidx(i), (b-1)*nso+cols+l, md+1) + sum(ff(1:p,e)*Lc(:,l))
              end do
            end do
          end do
        end do
      end do
    end do
    deallocate(tin, joa, ci, zc, zcn, nidx, znear, znearn, ikq, ikc, C1, C2, C3, C4, As, Ad, A1, A2, A3, A4, Gcq, Gpc, ff)
  end subroutine axissymstok_slpn_blockmat_nmode_r64

  subroutine axissymstok_dlp_blockmat_r64(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, mu, A)
    integer(8),   intent(in)    :: nt, p, np, iside, iclosed
    complex(r64), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: sws(p*np), tpan(np+1), mu
    real(r64),    intent(inout) :: A(2*nt, 2*np*p)
    integer(8) :: q, nso, i, j, jp, iq, ia, l, nc, cjo, k, jo, ne, npa, nq_az
    real(r64)  :: twopi, sumws, rho, rhop, zh, drho, rr2, chi, ws, spd, rt, zti, nrp, nzp, vk, ve, Fn, An, dFn
    real(r64)  :: t1, tN, tm, denom, rlo, rhi, tgi, split, thr, vrr, vrz, vzr, vzz
    complex(r64) :: ic, za, zb, Dz
    real(r64)    :: tglp(p), wglp(p), Dp(p,p), tglq(2*p), wglq(2*p), Dq(2*p,2*p), IP2(2*p,p)
    real(r64)    :: tc(p,np), rk(p), re2(2), IPk(p,p), IPe2(2,p), wsq(2*p), wsp(p)
    complex(r64) :: Yp(p), Yq(2*p), dYp(p), dYq(2*p), nvq(2*p), nvp(p), wxpq(2*p), ec2(2), xc(p,np)
    real(r64),    allocatable :: tin(:), tt_az(:), ww_az(:)
    logical,      allocatable :: cl(:)
    integer(8),   allocatable :: ci(:)
    complex(r64), allocatable :: zc(:)
    real(r64),    allocatable :: C1(:,:), C2(:,:), As(:,:), A1(:,:), A2(:,:), A3(:,:), A4(:,:), Gcl(:,:), Gp(:,:), ff(:,:)
    complex(r64), allocatable :: C3(:,:), C4(:,:), Ad(:,:)
    q = 2*p; nso = np*p; ic = (0.0_r64,1.0_r64); twopi = 2.0_r64*acos(-1.0_r64)
    call gauss_r64(p,     tglp, wglp, Dp)
    call gauss_r64(2_8*p, tglq, wglq, Dq)
    call lagrange_interp_r64(p, tglp, q, tglq, IP2)
    do j = 1, np
      do i = 1, p
        xc(i,j) = sx((j-1)*p+i)
        tc(i,j) = tpan(j) + (1.0_r64+tglp(i))/2.0_r64*(tpan(j+1)-tpan(j))
      end do
    end do
    t1 = tc(1,1); tN = tc(p,np)
    ! ---- dyadic pole refinement of tpan -> tin (mirror StoDLPAxiBlockMat_v0: 1st/last only) ----
    allocate(tin(np+1+256)); ne = np+1; tin(1:ne) = tpan(1:np+1)
    do while (tin(2) > t1)
      split = 0.5_r64*(tin(1)+tin(2))
      do i = ne, 2, -1; tin(i+1) = tin(i); end do
      tin(2) = split; ne = ne+1
    end do
    do while (tin(ne-1) < tN)
      split = 0.5_r64*(tin(ne)+tin(ne-1))
      tin(ne+1) = tin(ne); tin(ne) = split; ne = ne+1
    end do
    npa = ne-1
    thr = 0.1_r64
    call weight_setup_r64(p, 32_8*p, 1.0e-7_r64, nq_az, tt_az, ww_az)   ! graded azimuthal mesh (far close pairs)
    allocate(cl(nt), ci(nt), zc(nt), C1(2*nt,2*q), C2(2*nt,2*q), C3(2*nt,2*q), C4(2*nt,2*q))
    allocate(As(q,nt), Ad(q,nt), A1(q,nt), A2(q,nt), A3(q,nt), A4(q,nt), Gcl(2*nt,2*q), Gp(nt,p), ff(p,4))
    A = 0.0_r64
    do k = 1, npa
      ! ---- this refined panel: parent jo, coarse->refined interp IPk, positions/endpoints ----
      tm = 0.5_r64*(tin(k)+tin(k+1)); jo = 1
      do j = 1, np
        if (tm >= tpan(j) .and. tm < tpan(j+1)) jo = j
      end do
      denom = tpan(jo+1) - tpan(jo)
      do i = 1, p
        tgi = tin(k) + (1.0_r64+tglp(i))/2.0_r64*(tin(k+1)-tin(k))
        rk(i) = (2.0_r64*tgi - (tpan(jo)+tpan(jo+1)))/denom
      end do
      call lagrange_interp_r64(p, tglp, p, rk, IPk)
      Yp = matmul(IPk, xc(:,jo))
      rlo = (2.0_r64*tin(k)   - (tpan(jo)+tpan(jo+1)))/denom
      rhi = (2.0_r64*tin(k+1) - (tpan(jo)+tpan(jo+1)))/denom
      re2(1) = rlo; re2(2) = rhi
      call lagrange_interp_r64(p, tglp, 2_8, re2, IPe2)
      ec2 = matmul(IPe2, xc(:,jo)); za = ec2(1); zb = ec2(2)
      cjo = (jo-1)*p
      ! ---- refined-panel geometry (upsample + speed/source normal via D-matrix) ----
      Yq = matmul(IP2, Yp); dYq = matmul(Dq, Yq)
      do iq = 1, q
        spd = abs(dYq(iq)); nvq(iq) = -ic*dYq(iq)/spd; wsq(iq) = wglq(iq)*spd; wxpq(iq) = dYq(iq)*wglq(iq)
      end do
      dYp = matmul(Dp, Yp)
      do jp = 1, p
        spd = abs(dYp(jp)); wsp(jp) = wglp(jp)*spd; nvp(jp) = -ic*dYp(jp)/spd
      end do
      sumws = sum(wsp)
      nc = 0
      do i = 1, nt
        cl(i) = (abs(tx(i)-za) + abs(tx(i)-zb)) < 1.50_r64*sumws !1.85 causes early stagation for some reason... Hai
        if (cl(i)) then; nc = nc + 1; ci(nc) = i; zc(nc) = tx(i); end if
      end do
      ! ---- near: 2p DLP close (coef + sdspecialquad, src normal nvq), fold q->p (IP2) -> coarse (IPk) ----
      if (nc > 0) then
        call axissymstok_dlp_coef_r64(nc, zc(1:nc), q, Yq, nvq, mu, &
             C1(1:2*nc,1:2*q), C2(1:2*nc,1:2*q), C3(1:2*nc,1:2*q), C4(1:2*nc,1:2*q))
        call sdspecialquad_r64(nc, zc(1:nc), q, Yq, nvq, wxpq, za, zb, iside, &
             As(1:q,1:nc), Ad(1:q,1:nc), A1(1:q,1:nc), A2(1:q,1:nc), A3(1:q,1:nc), A4(1:q,1:nc))
        do iq = 1, q
          do ia = 1, nc
            Dz = cmplx(A1(iq,ia), -A2(iq,ia), r64)
            Gcl(ia,    iq)   = twopi*C1(ia,iq)*As(iq,ia)     + C2(ia,iq)*wsq(iq) &
                 + twopi*real(Ad(iq,ia)*(C3(ia,iq)/nvq(iq)), r64) - twopi*real(Dz*(C4(ia,iq)/nvq(iq)), r64)
            Gcl(ia,    q+iq) = twopi*C1(ia,q+iq)*As(iq,ia)   + C2(ia,q+iq)*wsq(iq) &
                 + twopi*real(Ad(iq,ia)*(C3(ia,q+iq)/nvq(iq)), r64) - twopi*real(Dz*(C4(ia,q+iq)/nvq(iq)), r64)
            Gcl(nc+ia, iq)   = twopi*C1(nc+ia,iq)*As(iq,ia)  + C2(nc+ia,iq)*wsq(iq) &
                 + twopi*real(Ad(iq,ia)*(C3(nc+ia,iq)/nvq(iq)), r64) - twopi*real(Dz*(C4(nc+ia,iq)/nvq(iq)), r64)
            Gcl(nc+ia, q+iq) = twopi*C1(nc+ia,q+iq)*As(iq,ia)+ C2(nc+ia,q+iq)*wsq(iq) &
                 + twopi*real(Ad(iq,ia)*(C3(nc+ia,q+iq)/nvq(iq)), r64) - twopi*real(Dz*(C4(nc+ia,q+iq)/nvq(iq)), r64)
          end do
        end do
        Gp(1:nc,:) = matmul(Gcl(1:nc, 1:q), IP2)
        do ia = 1, nc; do l = 1, p; A(ci(ia),    cjo+l)     = A(ci(ia),    cjo+l)     + sum(Gp(ia,1:p)*IPk(:,l)); end do; end do
        Gp(1:nc,:) = matmul(Gcl(1:nc, q+1:2*q), IP2)
        do ia = 1, nc; do l = 1, p; A(ci(ia),    nso+cjo+l) = A(ci(ia),    nso+cjo+l) + sum(Gp(ia,1:p)*IPk(:,l)); end do; end do
        Gp(1:nc,:) = matmul(Gcl(nc+1:2*nc, 1:q), IP2)
        do ia = 1, nc; do l = 1, p; A(nt+ci(ia), cjo+l)     = A(nt+ci(ia), cjo+l)     + sum(Gp(ia,1:p)*IPk(:,l)); end do; end do
        Gp(1:nc,:) = matmul(Gcl(nc+1:2*nc, q+1:2*q), IP2)
        do ia = 1, nc; do l = 1, p; A(nt+ci(ia), nso+cjo+l) = A(nt+ci(ia), nso+cjo+l) + sum(Gp(ia,1:p)*IPk(:,l)); end do; end do
      end if
      ! ---- far: dlp9 meridian (plus) on the refined p nodes, src normal nvp, fold via IPk ----
      do i = 1, nt
        if (cl(i)) cycle
        rt = real(tx(i),r64); zti = aimag(tx(i))
        do jp = 1, p
          rhop = real(Yp(jp),r64); zh = zti - aimag(Yp(jp)); nrp = real(nvp(jp),r64); nzp = aimag(nvp(jp))
          rho = rt; drho = rho - rhop; rr2 = drho*drho + zh*zh
          if (sqrt(rr2) < thr) then
            ! close pair in the far bucket (target near the source ring): cancellation-free azimuthal quad
            call axissymstok_dlp_aziquad_r64(rho, rhop, zh, nrp, nzp, nq_az, tt_az, ww_az, vrr, vrz, vzr, vzz)
            ff(jp,1) = vrr*rhop*wsp(jp); ff(jp,2) = vrz*rhop*wsp(jp)
            ff(jp,3) = vzr*rhop*wsp(jp); ff(jp,4) = vzz*rhop*wsp(jp)
            cycle
          end if
          chi = 1.0_r64 + rr2/(2.0_r64*rho*rhop)
          call carrier_r64(chi, 0_8, vk, ve, Fn, An, dFn)
          block
            real(r64) :: t2,t3,t4,t5,t6,t7,t8,t9,t11,t13,t15,t16,t17,t18,t19,t20,t21,t23,t24,t25,t26,t27
            real(r64) :: t30,t31,t32,t33,t39,t40,t41,t43,t44,t45,t46,t47,t48,t49,t50,t57,t58,t59,t62,t65,t70
            real(r64) :: c_1,c_2,c_11,c_12,c_31,c_32,c_41,c_42
            t2=nrp*rho; t3=nrp*rhop; t4=rho*rhop; t5=nzp*zh; t6=chi+1.0_r64; t7=chi*chi; t8=chi**3
            t11=rho*rho; t13=rhop*rhop; t15=zh*zh; t17=1.0_r64/acos(-1.0_r64); t20=chi-1.0_r64; t21=-rhop
            t9=t7*t7; t16=t2*3.0_r64; t18=rho*t5; t19=rhop*t5; t23=-t3; t24=-t5; t25=rho*t2; t26=rhop*t3
            t27=t3*t13; t30=chi*rhop*t2*2.0_r64; t31=t5*t13; t32=1.0_r64/t6
            t33=t7-1.0_r64; t39=rho*t16; t40=t11*t16; t41=t26*3.0_r64; t43=-rhop*t2*4.0_r64; t44=-t18; t45=-t19; t46=t4*t16
            t47=t2*t13*7.0_r64; t48=1.0_r64/t20; t49=t48*t48; t50=-t4*t5*4.0_r64; t57=t13*t24; t58=rhop*t2*t7*2.0_r64
            t59=-1.0_r64; t62=1.0_r64/t4**2.5_r64; t65=1.0_r64/t33; t70=-1.0_r64
            c_1  = rhop*t17*t62*t65*(t40+t47+t50+chi*t31-chi*t2*t4*11.0_r64+chi*t5*t11+chi*t13*t23 &
                 +t2*t4*t8*8.0_r64+t4*t5*t7*2.0_r64-t2*t7*t11*2.0_r64-t2*t7*t13*4.0_r64)*(-1.0_r64/8.0_r64)
            c_2  = (rhop*t17*t32*t49*t62*t70*(t27*3.0_r64-t31*3.0_r64+t46-t5*t11*3.0_r64+t7*t27+t7*t57 &
                 +chi*t4*t5*10.0_r64-chi*t2*t11*6.0_r64-chi*t2*t13*16.0_r64+t2*t4*t7*17.0_r64-t2*t4*t9*8.0_r64 &
                 -t4*t5*t8*2.0_r64+t2*t8*t11*2.0_r64+t2*t8*t13*4.0_r64+t7*t11*t24))/8.0_r64
            c_11 = rhop*t17*t62*t65*zh*(t18+t43+t58+chi*t25+chi*t26+chi*t45)*(-1.0_r64/8.0_r64)
            c_12 = (rhop*t17*t32*t49*t62*zh*(t19*(-3.0_r64)+t39+t41+chi*t18*4.0_r64+t7*t25+t7*t26+t7*t45 &
                 -chi*rhop*t2*10.0_r64+rhop*t2*t8*2.0_r64))/8.0_r64
            c_31 = (rhop*t17*t62*t65*zh*(t19-t25*3.0_r64+t30+chi*t44+t3*t21+t7*t25*2.0_r64))/8.0_r64
            c_32 = (rhop*t17*t32*t49*t62*zh*(t18*3.0_r64-t58-chi*t19*4.0_r64+chi*t25*6.0_r64+chi*t26*4.0_r64 &
                 -rhop*t2*6.0_r64+t7*t18-t8*t25*2.0_r64))/8.0_r64
            c_41 = (rhop*t15*t17*t59*t62*t65*(t5+t23+chi*t2))/8.0_r64
            c_42 = (rhop*t15*t17*t32*t49*t62*(t16-chi*t3*4.0_r64+chi*t5*4.0_r64+t2*t7))/8.0_r64
            ff(jp,1) = (c_1 *vk + c_2 *ve)*wsp(jp); ff(jp,2) = (c_11*vk + c_12*ve)*wsp(jp)
            ff(jp,3) = (c_31*vk + c_32*ve)*wsp(jp); ff(jp,4) = (c_41*vk + c_42*ve)*wsp(jp)
          end block
        end do
        do l = 1, p
          A(i,     cjo+l)     = A(i,     cjo+l)     + sum(ff(1:p,1)*IPk(:,l))
          A(i,     nso+cjo+l) = A(i,     nso+cjo+l) + sum(ff(1:p,2)*IPk(:,l))
          A(nt+i,  cjo+l)     = A(nt+i,  cjo+l)     + sum(ff(1:p,3)*IPk(:,l))
          A(nt+i,  nso+cjo+l) = A(nt+i,  nso+cjo+l) + sum(ff(1:p,4)*IPk(:,l))
        end do
      end do
    end do
    deallocate(tin, tt_az, ww_az, cl, ci, zc, C1, C2, C3, C4, As, Ad, A1, A2, A3, A4, Gcl, Gp, ff)
  end subroutine axissymstok_dlp_blockmat_r64

  subroutine axissymstok_dlp_blockmat_nmode_r64(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, mu, A)
    integer(8),   intent(in)    :: nt, p, np, M, iside, iclosed
    complex(r64), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: sws(p*np), tpan(np+1), mu
    complex(r64), intent(inout) :: A(3*nt, 3*np*p, M+1)
    integer(8), parameter :: Ksub = 4
    integer(8) :: q, pq, nso, i, j, jp, iq, ia, l, e, md, ar, b, c, qo, ne, npa, nk, nkc, jj, cols
    real(r64)  :: twopi, sumwso, sumwsc, rho, rhop, zh, rr2, chi, ws, spd, rt, zti, nrp, nzp, vk, ve, Fn, An, dFn, rn
    real(r64)  :: st1, stN, tm, denom, rlo, rhi, tgi, split
    complex(r64) :: ic, zac, zbc, cc1, cc2, cc3, cc4, Slog, Dval, Dz, Acau, Ahy
    real(r64)    :: tglp(p), wglp(p), Dp(p,p), tglq(2*p), wglq(2*p), Dq(2*p,2*p), IP2(2*p,p)
    real(r64)    :: tc(p,np), rk(p), re2(2), Lc(p,p), IPe2(2,p), IPqc(2*p,p), wsq(2*p), wsp(p)
    complex(r64) :: Ypc(p), Ypb(p), Yq(2*p), dYq(2*p), dYp(p), nvq(2*p), wxpq(2*p), ec2(2), xc(p,np)
    real(r64),    allocatable :: tin(:)
    integer(8),   allocatable :: joa(:), ci(:), nidx(:)
    logical,      allocatable :: ikq(:), ikc(:)
    complex(r64), allocatable :: zc(:), znear(:), C1(:,:,:), C2(:,:,:), C3(:,:,:), C4(:,:,:), Gcq(:,:), Gpc(:,:), ff(:,:)
    real(r64),    allocatable :: As(:,:), A1(:,:), A2(:,:), A3(:,:), A4(:,:)
    complex(r64), allocatable :: Ad(:,:)
    q = 2*p; nso = np*p; ic = (0.0_r64,1.0_r64); twopi = 2.0_r64*acos(-1.0_r64)
    call gauss_r64(p,     tglp, wglp, Dp)
    call gauss_r64(2_8*p, tglq, wglq, Dq)
    call lagrange_interp_r64(p, tglp, q, tglq, IP2)
    do j = 1, np
      do i = 1, p
        xc(i,j) = sx((j-1)*p+i)
        tc(i,j) = tpan(j) + (1.0_r64+tglp(i))/2.0_r64*(tpan(j+1)-tpan(j))
      end do
    end do
    st1 = tc(1,1); stN = tc(p,np)
    ! ---- two-level aux mesh: Ksub uniform subpanels per panel + dyadic pole refine ----
    allocate(tin(np*Ksub+1+512)); ne = 0
    do j = 1, np
      do i = 0, Ksub-1
        ne = ne+1; tin(ne) = tpan(j) + (tpan(j+1)-tpan(j))*real(i,r64)/real(Ksub,r64)
      end do
    end do
    ne = ne+1; tin(ne) = tpan(np+1)
    do while (tin(2) > st1)
      split = 0.5_r64*(tin(1)+tin(2)); do i = ne, 2, -1; tin(i+1) = tin(i); end do; tin(2) = split; ne = ne+1
    end do
    do while (tin(ne-1) < stN)
      split = 0.5_r64*(tin(ne)+tin(ne-1)); tin(ne+1) = tin(ne); tin(ne) = split; ne = ne+1
    end do
    npa = ne-1
    allocate(joa(npa))
    do c = 1, npa
      tm = 0.5_r64*(tin(c)+tin(c+1)); qo = 1
      do j = 1, np
        if (tm >= tpan(j) .and. tm < tpan(j+1)) qo = j
      end do
      joa(c) = qo
    end do
    allocate(ci(nt), zc(nt), nidx(nt), znear(nt), ikq(nt), ikc(nt))
    allocate(C1(3*nt,3*q,M+1), C2(3*nt,3*q,M+1), C3(3*nt,3*q,M+1), C4(3*nt,3*q,M+1))
    allocate(As(q,nt), Ad(q,nt), A1(q,nt), A2(q,nt), A3(q,nt), A4(q,nt), Gcq(nt,q), Gpc(nt,p), ff(p,9))
    A = (0.0_r64,0.0_r64)
    do pq = 1, np                                                   ! ===== outer level: original panel pq =====
      cols = (pq-1)*p
      do i = 1, p
        Ypc(i) = xc(i,pq)
      end do
      sumwso = 0.0_r64
      do jp = 1, p
        sumwso = sumwso + sws(cols+jp)
      end do
      nk = 0
      do i = 1, nt
        ikq(i) = (abs(tx(i)-sxlo(pq)) + abs(tx(i)-sxhi(pq))) < 1.85_r64*sumwso
        if (ikq(i)) then; nk = nk + 1; nidx(nk) = i; znear(nk) = tx(i); end if
      end do
      ! ---- outer far: naive dlp9far on the COARSE original nodes (provided snx, sws) ----
      do i = 1, nt
        if (ikq(i)) cycle
        rt = real(tx(i),r64); zti = aimag(tx(i))
        do jp = 1, p
          jj = cols + jp; rhop = real(sx(jj),r64); zh = zti - aimag(sx(jj)); ws = sws(jj)
          nrp = real(snx(jj),r64); nzp = aimag(snx(jj)); rho = rt
          rr2 = (rho-rhop)**2 + zh*zh; chi = 1.0_r64 + rr2/(2.0_r64*rho*rhop)
          do md = 0, M
            rn = real(md,r64); call carrier_r64(chi, md, vk, ve, Fn, An, dFn)
            block
              real(r64) :: t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t13,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26
              real(r64) :: t27,t28,t29,t30,t31,t32,t33,t34,t39,t40,t41,t43,t44,t45,t46,t47,t48,t49,t50,t57,t58,t59,t62,t65,t70
              complex(r64) :: t35,t36,t37,t38,t51,t52,t53,t54,t56,t60,t61,t67,t68,t69,t73, coK(9), coE(9)
              t2=nrp*rho; t3=nrp*rhop; t4=rho*rhop; t5=nzp*zh; t6=chi+1.0_r64; t7=chi*chi; t8=chi**3
              t10=rn*rn; t11=rho*rho; t13=rhop*rhop; t15=zh*zh; t17=1.0_r64/acos(-1.0_r64); t20=chi-1.0_r64; t21=-rhop
              t9=t7*t7; t16=t2*3.0_r64; t18=rho*t5; t19=rhop*t5; t22=t10*4.0_r64; t23=-t3; t24=-t5; t25=rho*t2; t26=rhop*t3
              t27=t3*t13; t28=chi*t2*4.0_r64; t29=rhop*t2*4.0_r64; t30=chi*rhop*t2*2.0_r64; t31=t5*t13; t32=1.0_r64/t6
              t33=t7-1.0_r64; t34=t4*t5*4.0_r64; t39=rho*t16; t40=t11*t16; t41=t26*3.0_r64; t43=-t29; t44=-t18; t45=-t19
              t46=t4*t16; t47=t2*t13*7.0_r64; t48=1.0_r64/t20; t49=t48*t48; t50=-t34; t57=t13*t24; t58=rhop*t2*t7*2.0_r64
              t59=t22-1.0_r64; t62=1.0_r64/t4**2.5_r64; t65=1.0_r64/t33; t70=1.0_r64/t59
              t35=t2*ic; t36=t3*ic; t37=t3*2.0_r64*ic; t38=t5*ic; t51=-t36; t52=chi*t35; t53=chi*t2*2.0_r64*ic
              t54=t19*ic; t56=rhop*t53; t60=t25*ic; t61=t26*ic; t67=rhop*t2*t10*4.0_r64*ic; t68=t10*t19*4.0_r64*ic
              t69=t10*t26*4.0_r64*ic; t73=t38+t51+t52
              coK(1) = cmplx(rhop*t17*t62*t65*(t40+t47+t50+chi*t31-chi*t2*t4*11.0_r64+chi*t5*t11+chi*t13*t23 &
                     -chi*t10*t31*4.0_r64+chi*t22*t27+t2*t4*t8*8.0_r64+t4*t5*t7*2.0_r64-t2*t7*t11*2.0_r64-t2*t7*t13*4.0_r64 &
                     -t2*t10*t13*4.0_r64+t4*t5*t22+chi*t2*t4*t10*8.0_r64-chi*t5*t10*t11*4.0_r64-t2*t7*t10*t11*4.0_r64 &
                     -t2*t7*t10*t13*8.0_r64+t2*t4*t8*t22+t4*t5*t7*t22)*(-1.0_r64/8.0_r64), 0.0_r64, r64)
              coE(1) = cmplx((rhop*t17*t32*t49*t62*t70*(t27*3.0_r64-t31*3.0_r64+t46-t5*t11*3.0_r64+t7*t27-t10*t27*12.0_r64 &
                     +t10*t31*12.0_r64+t7*t57+chi*t4*t5*10.0_r64-chi*t2*t11*6.0_r64-chi*t2*t13*16.0_r64+t2*t4*t7*17.0_r64 &
                     -t2*t4*t9*8.0_r64-t2*t4*t10*24.0_r64-t4*t5*t8*2.0_r64+t2*t8*t11*2.0_r64+t2*t8*t13*4.0_r64 &
                     +t5*t10*t11*12.0_r64+t7*t11*t24-t7*t10*t27*4.0_r64+t7*t22*t31-chi*t4*t5*t10*40.0_r64 &
                     +chi*t2*t10*t11*24.0_r64+chi*t2*t10*t13*64.0_r64-t2*t4*t7*t10*44.0_r64+t2*t4*t9*t10*20.0_r64 &
                     +t4*t5*t8*t10*8.0_r64-t2*t8*t10*t11*8.0_r64-t2*t8*t10*t13*16.0_r64+t5*t7*t11*t22))/8.0_r64, 0.0_r64, r64)
              coK(2) = rn/t4**1.5_r64*t17*(t54+t56-t60-t61)*(-3.0_r64/4.0_r64)
              coE(2) = (rn/t4**1.5_r64*t17*t48*t70*(t18*ic+t67-chi*t19*ic+chi*t60+chi*t61+chi*t68-rhop*t2*4.0_r64*ic &
                     -t10*t18*4.0_r64*ic+t7*t67-chi*t10*t25*4.0_r64*ic-chi*t10*t26*4.0_r64*ic+rhop*t2*t7*2.0_r64*ic))/4.0_r64
              coK(3) = cmplx(rhop*t17*t62*t65*zh*(t18+t43+t58+chi*t25+chi*t26+chi*t45-t10*t18*4.0_r64-chi*t10*t25*4.0_r64 &
                     -chi*t10*t26*4.0_r64+chi*t19*t22+rhop*t2*t22+rhop*t2*t7*t22)*(-1.0_r64/8.0_r64), 0.0_r64, r64)
              coE(3) = cmplx((rhop*t17*t32*t49*t62*zh*(t19*(-3.0_r64)+t39+t41+chi*t18*4.0_r64+t7*t25+t7*t26+t7*t45 &
                     -chi*rhop*t2*10.0_r64+rhop*t2*t8*2.0_r64))/8.0_r64, 0.0_r64, r64)
              coK(4) = rn*rhop/t4**1.5_r64*t17*(-t37+t38+t53)*(3.0_r64/4.0_r64)
              coE(4) = rn*t13*t17*t48*t62*t70*(t25*(-3.0_r64)*ic+t54+t56-t61-t68+t69-chi*t18*ic+t7*t25*2.0_r64*ic &
                     +chi*t10*t18*4.0_r64*ic+t7*t10*t25*4.0_r64*ic-chi*rhop*t2*t10*8.0_r64*ic)*(-1.0_r64/4.0_r64)
              coK(5) = cmplx(rhop/t4**1.5_r64*t17*(t5+t23+t28-t3*t10*2.0_r64+t5*t10*2.0_r64+chi*t2*t10*2.0_r64)*(-1.0_r64/4.0_r64), 0.0_r64, r64)
              coE(5) = cmplx((rhop/t4**1.5_r64*t17*t48*t70*(t16+chi*t3+chi*t24-t2*t7*4.0_r64-t2*t10*6.0_r64 &
                     -chi*t3*t10*4.0_r64+chi*t5*t22+t2*t7*t10*10.0_r64))/4.0_r64, 0.0_r64, r64)
              coK(6) = rn*t2*t13*t17*t62*zh*cmplx(0.0_r64,0.75_r64,r64)
              coE(6) = rn*t13*t17*t48*t62*t73*zh*(-1.0_r64/4.0_r64)
              coK(7) = cmplx((rhop*t17*t62*t65*zh*(t19-t25*3.0_r64+t30+chi*t44+t3*t21-t10*t19*4.0_r64+t7*t25*2.0_r64 &
                     +t22*t26+chi*t18*t22+t7*t22*t25-chi*rhop*t2*t10*8.0_r64))/8.0_r64, 0.0_r64, r64)
              coE(7) = cmplx((rhop*t17*t32*t49*t62*zh*(t18*3.0_r64-t58-chi*t19*4.0_r64+chi*t25*6.0_r64+chi*t26*4.0_r64 &
                     -rhop*t2*6.0_r64+t7*t18-t8*t25*2.0_r64))/8.0_r64, 0.0_r64, r64)
              coK(8) = rn*t2/t4**1.5_r64*t17*zh*cmplx(0.0_r64,0.75_r64,r64)
              coE(8) = rn/t4**1.5_r64*t17*t48*t73*zh*(-1.0_r64/4.0_r64)
              coK(9) = cmplx((rhop*t15*t17*t59*t62*t65*(t5+t23+chi*t2))/8.0_r64, 0.0_r64, r64)
              coE(9) = cmplx((rhop*t15*t17*t32*t49*t62*(t16-chi*t3*4.0_r64+chi*t5*4.0_r64+t2*t7))/8.0_r64, 0.0_r64, r64)
              do e = 1, 9
                ar = (e-1)/3 + 1; b = mod(e-1,3) + 1
                A((ar-1)*nt+i, (b-1)*nso+jj, md+1) = (coK(e)*vk + coE(e)*ve)*ws
              end do
            end block
          end do
        end do
      end do
      if (nk == 0) cycle
      ! ---- outer near: recurse into the Ksub(+dyadic) subpanels of q, inner re-detect ----
      do c = 1, npa
        if (joa(c) /= pq) cycle
        denom = tpan(pq+1) - tpan(pq)
        do i = 1, p
          tgi = tin(c) + (1.0_r64+tglp(i))/2.0_r64*(tin(c+1)-tin(c))
          rk(i) = (2.0_r64*tgi - (tpan(pq)+tpan(pq+1)))/denom
        end do
        call lagrange_interp_r64(p, tglp, p, rk, Lc)
        Ypb = matmul(Lc, xc(:,pq))
        rlo = (2.0_r64*tin(c)   - (tpan(pq)+tpan(pq+1)))/denom
        rhi = (2.0_r64*tin(c+1) - (tpan(pq)+tpan(pq+1)))/denom
        re2(1) = rlo; re2(2) = rhi
        call lagrange_interp_r64(p, tglp, 2_8, re2, IPe2)
        ec2 = matmul(IPe2, xc(:,pq)); zac = ec2(1); zbc = ec2(2)
        IPqc = matmul(IP2, Lc)
        Yq = matmul(IP2, Ypb); dYq = matmul(Dq, Yq)
        do iq = 1, q
          spd = abs(dYq(iq)); nvq(iq) = -ic*dYq(iq)/spd; wsq(iq) = wglq(iq)*spd; wxpq(iq) = dYq(iq)*wglq(iq)
        end do
        dYp = matmul(Dp, Ypb)
        do jp = 1, p
          wsp(jp) = wglp(jp)*abs(dYp(jp))
        end do
        sumwsc = sum(wsp)
        nkc = 0
        do i = 1, nk
          ikc(i) = (abs(znear(i)-zac) + abs(znear(i)-zbc)) < 1.5_r64*sumwsc
          if (ikc(i)) then; nkc = nkc + 1; ci(nkc) = i; zc(nkc) = znear(i); end if
        end do
        ! ---- inner near: 2p DLP close (coef + sdspecialquad), fold q->p_sub (IP2) -> p_coarse (Lc) ----
        if (nkc > 0) then
          call axissymstok_dlp_coef_nmode_r64(nkc, zc(1:nkc), q, Yq, nvq, M, mu, &
               C1(1:3*nkc,1:3*q,:), C2(1:3*nkc,1:3*q,:), C3(1:3*nkc,1:3*q,:), C4(1:3*nkc,1:3*q,:))
          call sdspecialquad_r64(nkc, zc(1:nkc), q, Yq, nvq, wxpq, zac, zbc, iside, &
               As(1:q,1:nkc), Ad(1:q,1:nkc), A1(1:q,1:nkc), A2(1:q,1:nkc), A3(1:q,1:nkc), A4(1:q,1:nkc))
          do md = 0, M
            do e = 1, 9
              ar = (e-1)/3 + 1; b = mod(e-1,3) + 1
              do iq = 1, q
                do ia = 1, nkc
                  cc1 = C1(3*(ia-1)+ar, 3*(iq-1)+b, md+1); cc2 = C2(3*(ia-1)+ar, 3*(iq-1)+b, md+1)
                  cc3 = C3(3*(ia-1)+ar, 3*(iq-1)+b, md+1); cc4 = C4(3*(ia-1)+ar, 3*(iq-1)+b, md+1)
                  Slog = cmplx(As(iq,ia),0.0_r64,r64); Dval = Ad(iq,ia); Dz = cmplx(A1(iq,ia),-A2(iq,ia),r64)
                  if (mod(e,2) == 0) then
                    Acau =  ic*twopi*aimag(Dval*(cc3/nvq(iq)))
                    Ahy  = -ic*twopi*aimag(Dz  *(cc4/nvq(iq)))
                  else
                    Acau = cmplx(twopi*real(Dval*(cc3/nvq(iq)), r64), 0.0_r64, r64)
                    Ahy  = cmplx(-twopi*real(Dz *(cc4/nvq(iq)), r64), 0.0_r64, r64)
                  end if
                  Gcq(ia,iq) = twopi*cc1*Slog + cc2*wsq(iq) + Acau + Ahy
                end do
              end do
              Gpc(1:nkc,1:p) = matmul(Gcq(1:nkc,1:q), IPqc)
              do ia = 1, nkc
                do l = 1, p
                  A((ar-1)*nt+nidx(ci(ia)), (b-1)*nso+cols+l, md+1) = &
                    A((ar-1)*nt+nidx(ci(ia)), (b-1)*nso+cols+l, md+1) + Gpc(ia,l)
                end do
              end do
            end do
          end do
        end if
        ! ---- inner far: naive dlp9far on the FINER subpanel nodes, fold via Lc ----
        do i = 1, nk
          if (ikc(i)) cycle
          rt = real(znear(i),r64); zti = aimag(znear(i))
          do md = 0, M
            rn = real(md,r64)
            do jp = 1, p
              rhop = real(Ypb(jp),r64); zh = zti - aimag(Ypb(jp)); rho = rt
              nrp = real(-ic*dYp(jp)/abs(dYp(jp)),r64); nzp = aimag(-ic*dYp(jp)/abs(dYp(jp)))
              rr2 = (rho-rhop)**2 + zh*zh; chi = 1.0_r64 + rr2/(2.0_r64*rho*rhop)
              call carrier_r64(chi, md, vk, ve, Fn, An, dFn)
              block
                real(r64) :: t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t13,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26
                real(r64) :: t27,t28,t29,t30,t31,t32,t33,t34,t39,t40,t41,t43,t44,t45,t46,t47,t48,t49,t50,t57,t58,t59,t62,t65,t70
                complex(r64) :: t35,t36,t37,t38,t51,t52,t53,t54,t56,t60,t61,t67,t68,t69,t73, coK(9), coE(9)
                t2=nrp*rho; t3=nrp*rhop; t4=rho*rhop; t5=nzp*zh; t6=chi+1.0_r64; t7=chi*chi; t8=chi**3
                t10=rn*rn; t11=rho*rho; t13=rhop*rhop; t15=zh*zh; t17=1.0_r64/acos(-1.0_r64); t20=chi-1.0_r64; t21=-rhop
                t9=t7*t7; t16=t2*3.0_r64; t18=rho*t5; t19=rhop*t5; t22=t10*4.0_r64; t23=-t3; t24=-t5; t25=rho*t2; t26=rhop*t3
                t27=t3*t13; t28=chi*t2*4.0_r64; t29=rhop*t2*4.0_r64; t30=chi*rhop*t2*2.0_r64; t31=t5*t13; t32=1.0_r64/t6
                t33=t7-1.0_r64; t34=t4*t5*4.0_r64; t39=rho*t16; t40=t11*t16; t41=t26*3.0_r64; t43=-t29; t44=-t18; t45=-t19
                t46=t4*t16; t47=t2*t13*7.0_r64; t48=1.0_r64/t20; t49=t48*t48; t50=-t34; t57=t13*t24; t58=rhop*t2*t7*2.0_r64
                t59=t22-1.0_r64; t62=1.0_r64/t4**2.5_r64; t65=1.0_r64/t33; t70=1.0_r64/t59
                t35=t2*ic; t36=t3*ic; t37=t3*2.0_r64*ic; t38=t5*ic; t51=-t36; t52=chi*t35; t53=chi*t2*2.0_r64*ic
                t54=t19*ic; t56=rhop*t53; t60=t25*ic; t61=t26*ic; t67=rhop*t2*t10*4.0_r64*ic; t68=t10*t19*4.0_r64*ic
                t69=t10*t26*4.0_r64*ic; t73=t38+t51+t52
                coK(1) = cmplx(rhop*t17*t62*t65*(t40+t47+t50+chi*t31-chi*t2*t4*11.0_r64+chi*t5*t11+chi*t13*t23 &
                       -chi*t10*t31*4.0_r64+chi*t22*t27+t2*t4*t8*8.0_r64+t4*t5*t7*2.0_r64-t2*t7*t11*2.0_r64-t2*t7*t13*4.0_r64 &
                       -t2*t10*t13*4.0_r64+t4*t5*t22+chi*t2*t4*t10*8.0_r64-chi*t5*t10*t11*4.0_r64-t2*t7*t10*t11*4.0_r64 &
                       -t2*t7*t10*t13*8.0_r64+t2*t4*t8*t22+t4*t5*t7*t22)*(-1.0_r64/8.0_r64), 0.0_r64, r64)
                coE(1) = cmplx((rhop*t17*t32*t49*t62*t70*(t27*3.0_r64-t31*3.0_r64+t46-t5*t11*3.0_r64+t7*t27-t10*t27*12.0_r64 &
                       +t10*t31*12.0_r64+t7*t57+chi*t4*t5*10.0_r64-chi*t2*t11*6.0_r64-chi*t2*t13*16.0_r64+t2*t4*t7*17.0_r64 &
                       -t2*t4*t9*8.0_r64-t2*t4*t10*24.0_r64-t4*t5*t8*2.0_r64+t2*t8*t11*2.0_r64+t2*t8*t13*4.0_r64 &
                       +t5*t10*t11*12.0_r64+t7*t11*t24-t7*t10*t27*4.0_r64+t7*t22*t31-chi*t4*t5*t10*40.0_r64 &
                       +chi*t2*t10*t11*24.0_r64+chi*t2*t10*t13*64.0_r64-t2*t4*t7*t10*44.0_r64+t2*t4*t9*t10*20.0_r64 &
                       +t4*t5*t8*t10*8.0_r64-t2*t8*t10*t11*8.0_r64-t2*t8*t10*t13*16.0_r64+t5*t7*t11*t22))/8.0_r64, 0.0_r64, r64)
                coK(2) = rn/t4**1.5_r64*t17*(t54+t56-t60-t61)*(-3.0_r64/4.0_r64)
                coE(2) = (rn/t4**1.5_r64*t17*t48*t70*(t18*ic+t67-chi*t19*ic+chi*t60+chi*t61+chi*t68-rhop*t2*4.0_r64*ic &
                       -t10*t18*4.0_r64*ic+t7*t67-chi*t10*t25*4.0_r64*ic-chi*t10*t26*4.0_r64*ic+rhop*t2*t7*2.0_r64*ic))/4.0_r64
                coK(3) = cmplx(rhop*t17*t62*t65*zh*(t18+t43+t58+chi*t25+chi*t26+chi*t45-t10*t18*4.0_r64-chi*t10*t25*4.0_r64 &
                       -chi*t10*t26*4.0_r64+chi*t19*t22+rhop*t2*t22+rhop*t2*t7*t22)*(-1.0_r64/8.0_r64), 0.0_r64, r64)
                coE(3) = cmplx((rhop*t17*t32*t49*t62*zh*(t19*(-3.0_r64)+t39+t41+chi*t18*4.0_r64+t7*t25+t7*t26+t7*t45 &
                       -chi*rhop*t2*10.0_r64+rhop*t2*t8*2.0_r64))/8.0_r64, 0.0_r64, r64)
                coK(4) = rn*rhop/t4**1.5_r64*t17*(-t37+t38+t53)*(3.0_r64/4.0_r64)
                coE(4) = rn*t13*t17*t48*t62*t70*(t25*(-3.0_r64)*ic+t54+t56-t61-t68+t69-chi*t18*ic+t7*t25*2.0_r64*ic &
                       +chi*t10*t18*4.0_r64*ic+t7*t10*t25*4.0_r64*ic-chi*rhop*t2*t10*8.0_r64*ic)*(-1.0_r64/4.0_r64)
                coK(5) = cmplx(rhop/t4**1.5_r64*t17*(t5+t23+t28-t3*t10*2.0_r64+t5*t10*2.0_r64+chi*t2*t10*2.0_r64)*(-1.0_r64/4.0_r64), 0.0_r64, r64)
                coE(5) = cmplx((rhop/t4**1.5_r64*t17*t48*t70*(t16+chi*t3+chi*t24-t2*t7*4.0_r64-t2*t10*6.0_r64 &
                       -chi*t3*t10*4.0_r64+chi*t5*t22+t2*t7*t10*10.0_r64))/4.0_r64, 0.0_r64, r64)
                coK(6) = rn*t2*t13*t17*t62*zh*cmplx(0.0_r64,0.75_r64,r64)
                coE(6) = rn*t13*t17*t48*t62*t73*zh*(-1.0_r64/4.0_r64)
                coK(7) = cmplx((rhop*t17*t62*t65*zh*(t19-t25*3.0_r64+t30+chi*t44+t3*t21-t10*t19*4.0_r64+t7*t25*2.0_r64 &
                       +t22*t26+chi*t18*t22+t7*t22*t25-chi*rhop*t2*t10*8.0_r64))/8.0_r64, 0.0_r64, r64)
                coE(7) = cmplx((rhop*t17*t32*t49*t62*zh*(t18*3.0_r64-t58-chi*t19*4.0_r64+chi*t25*6.0_r64+chi*t26*4.0_r64 &
                       -rhop*t2*6.0_r64+t7*t18-t8*t25*2.0_r64))/8.0_r64, 0.0_r64, r64)
                coK(8) = rn*t2/t4**1.5_r64*t17*zh*cmplx(0.0_r64,0.75_r64,r64)
                coE(8) = rn/t4**1.5_r64*t17*t48*t73*zh*(-1.0_r64/4.0_r64)
                coK(9) = cmplx((rhop*t15*t17*t59*t62*t65*(t5+t23+chi*t2))/8.0_r64, 0.0_r64, r64)
                coE(9) = cmplx((rhop*t15*t17*t32*t49*t62*(t16-chi*t3*4.0_r64+chi*t5*4.0_r64+t2*t7))/8.0_r64, 0.0_r64, r64)
                do e = 1, 9
                  ff(jp,e) = (coK(e)*vk + coE(e)*ve)*wsp(jp)
                end do
              end block
            end do
            do e = 1, 9
              ar = (e-1)/3 + 1; b = mod(e-1,3) + 1
              do l = 1, p
                A((ar-1)*nt+nidx(i), (b-1)*nso+cols+l, md+1) = &
                  A((ar-1)*nt+nidx(i), (b-1)*nso+cols+l, md+1) + sum(ff(1:p,e)*Lc(:,l))
              end do
            end do
          end do
        end do
      end do
    end do
    deallocate(tin, joa, ci, zc, nidx, znear, ikq, ikc, C1, C2, C3, C4, As, Ad, A1, A2, A3, A4, Gcq, Gpc, ff)
  end subroutine axissymstok_dlp_blockmat_nmode_r64

  subroutine axissymstok_dlpn_blockmat_r64(nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, mu, A)
    integer(8),   intent(in)    :: nt, p, np, iside, iclosed
    complex(r64), intent(in)    :: tx(nt), tnx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: sws(p*np), tpan(np+1), mu
    real(r64),    intent(inout) :: A(2*nt, 2*np*p)
    integer(8) :: q, pnp, pnsa, i, j, jp, k, ne, npa, nk, nf, ia, ib, jo, l, nsc, nq_az
    real(r64)  :: sumws, dmin, dd, spd, t1, tN, tm, denom, tgi, split
    real(r64)  :: twopi, thr, drho, rr2, zti, rho, rhop, zh, nr, nz, nps, nzs, sc, Irr, Irz, Izr, Izz
    complex(r64) :: ic, zlo, zhi, Dz, Dzz
    real(r64)    :: tglp(p), wglp(p), Dp(p,p), tglq(2*p), wglq(2*p), Dq(2*p,2*p), IP2(2*p,p)
    real(r64)    :: rk(p), rk2(2), IPe2(2,p), IPk(p,p)
    complex(r64) :: xc(p,np), Yp(p), dYp(p), ec2(2)
    complex(r64) :: skx(p), sknx(p), skwxp(p), skfx(2*p), dY(2*p), csx(2*p), csnx(2*p), cswxp(2*p)
    real(r64)    :: skws(p), csws(2*p)
    real(r64),    allocatable :: tin(:), IM(:,:), Gc(:,:), Gf(:,:), Gfold(:,:), B(:,:), tt_az(:), ww_az(:)
    real(r64),    allocatable :: cc1(:,:), cc2(:,:), asq(:,:), a1q(:,:), a2q(:,:), a3q(:,:), a4q(:,:)
    complex(r64), allocatable :: cc3(:,:), cc4(:,:), cc5(:,:), adq(:,:)
    integer(8),   allocatable :: joa(:), ic_(:), ifr(:)
    logical,      allocatable :: ik(:)
    complex(r64), allocatable :: sax(:), sanx(:), sawxp(:), saxlo(:), saxhi(:)
    real(r64),    allocatable :: saws(:)
    complex(r64), allocatable :: tkx(:), tknx(:), tfx(:), tfnx(:)
    q = 2*p; pnp = p*np; ic = (0.0_r64,1.0_r64); twopi = 2.0_r64*acos(-1.0_r64); thr = 0.1_r64
    call weight_setup_r64(p, 32_8*p, 1.0e-7_r64, nq_az, tt_az, ww_az)   ! graded azimuthal mesh (close-in-far pairs)
    call gauss_r64(p,     tglp, wglp, Dp)
    call gauss_r64(2_8*p, tglq, wglq, Dq)
    call lagrange_interp_r64(p, tglp, q, tglq, IP2)
    do j = 1, np
      do i = 1, p
        xc(i,j) = sx((j-1)*p+i)
      end do
    end do
    t1 = tpan(1) + (1.0_r64+tglp(1))/2.0_r64*(tpan(2)-tpan(1))
    tN = tpan(np) + (1.0_r64+tglp(p))/2.0_r64*(tpan(np+1)-tpan(np))
    allocate(tin(np+1+256)); ne = np+1; tin(1:ne) = tpan(1:np+1)
    do while (tin(2) > t1)
      split = 0.5_r64*(tin(1)+tin(2))
      do i = ne, 2, -1; tin(i+1) = tin(i); end do
      tin(2) = split; ne = ne+1
    end do
    do while (tin(ne-1) < tN)
      split = 0.5_r64*(tin(ne)+tin(ne-1))
      tin(ne+1) = tin(ne); tin(ne) = split; ne = ne+1
    end do
    npa = ne-1; pnsa = p*npa
    allocate(joa(npa), IM(pnsa,pnp), sax(pnsa), sanx(pnsa), saws(pnsa), sawxp(pnsa), saxlo(npa), saxhi(npa))
    IM = 0.0_r64
    do k = 1, npa
      tm = 0.5_r64*(tin(k)+tin(k+1)); jo = 1
      do j = 1, np
        if (tm >= tpan(j) .and. tm < tpan(j+1)) jo = j
      end do
      joa(k) = jo; denom = tpan(jo+1) - tpan(jo)
      do i = 1, p
        tgi = tin(k) + (1.0_r64+tglp(i))/2.0_r64*(tin(k+1)-tin(k))
        rk(i) = (2.0_r64*tgi - (tpan(jo)+tpan(jo+1)))/denom
      end do
      call lagrange_interp_r64(p, tglp, p, rk, IPk)
      Yp = matmul(IPk, xc(:,jo)); dYp = matmul(Dp, Yp)
      do jp = 1, p
        spd = abs(dYp(jp))
        sax((k-1)*p+jp)   = Yp(jp)
        sanx((k-1)*p+jp)  = -ic*dYp(jp)/spd
        saws((k-1)*p+jp)  = wglp(jp)*spd
        sawxp((k-1)*p+jp) = dYp(jp)*wglp(jp)
        do l = 1, p; IM((k-1)*p+jp, (jo-1)*p+l) = IPk(jp,l); end do
      end do
      rk2(1) = (2.0_r64*tin(k)   - (tpan(jo)+tpan(jo+1)))/denom
      rk2(2) = (2.0_r64*tin(k+1) - (tpan(jo)+tpan(jo+1)))/denom
      call lagrange_interp_r64(p, tglp, 2_8, rk2, IPe2)
      ec2 = matmul(IPe2, xc(:,jo)); saxlo(k) = ec2(1); saxhi(k) = ec2(2)
    end do
    allocate(ik(nt), ic_(nt), ifr(nt), tkx(nt), tknx(nt), tfx(nt), tfnx(nt))
    allocate(Gc(2*nt,2*q), Gf(2*nt,2*p), Gfold(2*nt,2*p), B(2*nt,2*pnsa))
    allocate(cc1(2*nt,2*q), cc2(2*nt,2*q), cc3(2*nt,2*q), cc4(2*nt,2*q), cc5(2*nt,2*q))
    allocate(asq(q,nt), adq(q,nt), a1q(q,nt), a2q(q,nt), a3q(q,nt), a4q(q,nt))
    B = 0.0_r64
    do k = 1, npa
      do jp = 1, p
        skx(jp)  = sax((k-1)*p+jp);  sknx(jp)  = sanx((k-1)*p+jp)
        skws(jp) = saws((k-1)*p+jp); skwxp(jp) = sawxp((k-1)*p+jp)
      end do
      zlo = saxlo(k); zhi = saxhi(k); sumws = sum(skws)
      nk = 0; nf = 0
      do i = 1, nt
        if ((abs(tx(i)-zlo)+abs(tx(i)-zhi)) < 1.85_r64*sumws) then
          ik(i) = .true.;  nk = nk+1; ic_(nk) = i; tkx(nk) = tx(i); tknx(nk) = tnx(i)
        else
          ik(i) = .false.; nf = nf+1; ifr(nf) = i; tfx(nf) = tx(i); tfnx(nf) = tnx(i)
        end if
      end do
      if (nk > 0) then
        dmin = huge(1.0_r64)
        do ia = 1, nk; do jp = 1, p; dd = abs(tkx(ia)-skx(jp)); if (dd < dmin) dmin = dd; end do; end do
        if (dmin < 1.0e-3_r64*sumws) then
          nsc = q; skfx = matmul(IP2, skx); dY = matmul(Dq, skfx)
          do jp = 1, q
            spd = abs(dY(jp))
            csx(jp) = skfx(jp); csnx(jp) = -ic*dY(jp)/spd; csws(jp) = wglq(jp)*spd; cswxp(jp) = dY(jp)*wglq(jp)
          end do
        else
          nsc = p
          do jp = 1, p
            csx(jp) = skx(jp); csnx(jp) = sknx(jp); csws(jp) = skws(jp); cswxp(jp) = skwxp(jp)
          end do
        end if
        ! ---- near: coef + special quad + 5-bucket D' combine on the close (nsc) nodes ----
        call axissymstok_dlpn_coef_r64(nk, tkx, tknx, nsc, csx(1:nsc), csnx(1:nsc), mu, &
             cc1(1:2*nk,1:2*nsc), cc2(1:2*nk,1:2*nsc), cc3(1:2*nk,1:2*nsc), cc4(1:2*nk,1:2*nsc), cc5(1:2*nk,1:2*nsc))
        call sdspecialquad_r64(nk, tkx, nsc, csx(1:nsc), csnx(1:nsc), cswxp(1:nsc), zlo, zhi, iside, &
             asq(1:nsc,1:nk), adq(1:nsc,1:nk), a1q(1:nsc,1:nk), a2q(1:nsc,1:nk), a3q(1:nsc,1:nk), a4q(1:nsc,1:nk))
        do j = 1, nsc
          do i = 1, nk
            Dz = cmplx(a1q(j,i), -a2q(j,i), r64); Dzz = cmplx(a3q(j,i), -a4q(j,i), r64)
            Gc(i,    j)     = twopi*cc1(i,j)*asq(j,i)         + cc2(i,j)*csws(j) &
                 + twopi*real(adq(j,i)*(cc3(i,j)/csnx(j)),r64) - twopi*real(Dz*(cc4(i,j)/csnx(j)),r64) + twopi*real(Dzz*(cc5(i,j)/csnx(j)),r64)
            Gc(i,    nsc+j) = twopi*cc1(i,nsc+j)*asq(j,i)     + cc2(i,nsc+j)*csws(j) &
                 + twopi*real(adq(j,i)*(cc3(i,nsc+j)/csnx(j)),r64) - twopi*real(Dz*(cc4(i,nsc+j)/csnx(j)),r64) + twopi*real(Dzz*(cc5(i,nsc+j)/csnx(j)),r64)
            Gc(nk+i, j)     = twopi*cc1(nk+i,j)*asq(j,i)      + cc2(nk+i,j)*csws(j) &
                 + twopi*real(adq(j,i)*(cc3(nk+i,j)/csnx(j)),r64) - twopi*real(Dz*(cc4(nk+i,j)/csnx(j)),r64) + twopi*real(Dzz*(cc5(nk+i,j)/csnx(j)),r64)
            Gc(nk+i, nsc+j) = twopi*cc1(nk+i,nsc+j)*asq(j,i)  + cc2(nk+i,nsc+j)*csws(j) &
                 + twopi*real(adq(j,i)*(cc3(nk+i,nsc+j)/csnx(j)),r64) - twopi*real(Dz*(cc4(nk+i,nsc+j)/csnx(j)),r64) + twopi*real(Dzz*(cc5(nk+i,nsc+j)/csnx(j)),r64)
          end do
        end do
        if (nsc == q) then
          do ia = 1, 2*nk
            do jp = 1, p
              Gfold(ia,jp)   = sum(Gc(ia,1:q)*IP2(:,jp))
              Gfold(ia,p+jp) = sum(Gc(ia,q+1:2*q)*IP2(:,jp))
            end do
          end do
        else
          do ia = 1, 2*nk
            do jp = 1, p
              Gfold(ia,jp)   = Gc(ia,jp)
              Gfold(ia,p+jp) = Gc(ia,p+jp)
            end do
          end do
        end if
        do ia = 1, nk
          do jp = 1, p
            B(ic_(ia),    (k-1)*p+jp)      = Gfold(ia,    jp)
            B(ic_(ia),    pnsa+(k-1)*p+jp) = Gfold(ia,    p+jp)
            B(nt+ic_(ia), (k-1)*p+jp)      = Gfold(nk+ia, jp)
            B(nt+ic_(ia), pnsa+(k-1)*p+jp) = Gfold(nk+ia, p+jp)
          end do
        end do
      end if
      if (nf > 0) then
        ! ---- far: close-in-far pairs (dist<thr) use the cancellation-free graded azimuthal quad,
        !      the rest the analytic (K,E carrier) far; fold via IM ----
        do i = 1, nf
          rho = real(tfx(i),r64); zti = aimag(tfx(i)); nr = real(tfnx(i),r64); nz = aimag(tfnx(i))
          do j = 1, p
            rhop = real(skx(j),r64); zh = zti - aimag(skx(j)); nps = real(sknx(j),r64); nzs = aimag(sknx(j))
            drho = rho - rhop; rr2 = drho*drho + zh*zh
            if (sqrt(rr2) < thr) then
              call axissymstok_dlpn_aziquad_r64(rho, rhop, zh, nr, nz, nps, nzs, mu, nq_az, tt_az, ww_az, Irr, Irz, Izr, Izz)
            else
              block
                real(r64) :: chi, vk, ve, Fn, An, dFn
                real(r64) :: t2, t3, t4, t5, t6, t7, t8, t9, t10, t12, t14, t15
                real(r64) :: t16, t17, t18, t19, t20, t21, t22, t27, t28, t29, t30, t31
                real(r64) :: t32, t33, t37, t38, t39, t40, t41, t42, t45, t46, t48, t52
                real(r64) :: t53, t11, t13, t24, t26, t34, t35, t36, t47, t49, t50, t54
                real(r64) :: t55, t57, t59, t60, t61, t64, t65, t66, t72, t73, t74, t75
                real(r64) :: t76, t79, t81, t83, t84, t86, t87, t88, t93, t96, t97, t98
                real(r64) :: t103, t51, t56, t62, t67, t68, t69, t70, t77, t78, t80, t85
                real(r64) :: t89, t91, t92, t94, t95, t99, t102, t107, t114, t71, t90, t100
                real(r64) :: t101, t104, t109, t110, t113, t115, t121, t105, t108, t111, t116, t117
                real(r64) :: t118, t106, t112, t119, t120, t122, t123, et1, et2, et3, et4
                chi = (rho*rho+rhop*rhop+zh*zh)/(2.0_r64*rho*rhop)
                call carrier_r64(chi, 0_8, vk, ve, Fn, An, dFn)
                t2 = nr*rho
                t3 = nps*rhop
                t4 = rho*rhop
                t5 = nz*zh
                t6 = nzs*zh
                t7 = chi*2.0_r64
                t8 = chi+1.0_r64
                t9 = chi**2
                t10 = chi**3
                t12 = chi**5
                t14 = rho**2
                t15 = rho**3
                t16 = rhop**2
                t17 = rhop**3
                t18 = ve*2.0_r64
                t19 = ve*6.0_r64
                t20 = vk*2.0_r64
                t21 = vk*6.0_r64
                t22 = zh**2
                t27 = 1.0_r64/acos(-1.0_r64)
                t28 = chi-1.0_r64
                t29 = ve*1.8e+1_r64
                t30 = ve*3.0e+1_r64
                t31 = -vk
                t32 = vk*1.0e+1_r64
                t33 = vk*1.5e+1_r64
                t37 = chi*ve*8.0_r64
                t38 = chi*ve*1.2e+1_r64
                t39 = chi*ve*4.5e+1_r64
                t40 = chi*ve*5.8e+1_r64
                t41 = chi*vk*-2.0_r64
                t42 = chi*vk*-6.0_r64
                t45 = chi*vk*1.6e+1_r64
                t46 = chi*vk*2.0e+1_r64
                t48 = sqrt(2.0_r64)
                t52 = chi*vk*-1.0e+1_r64
                t53 = chi*vk*-1.5e+1_r64
                t11 = t9**2
                t13 = t9**3
                t24 = t7*vk
                t26 = t4*2.0_r64
                t34 = t7+2.0_r64
                t35 = nps*rho*t5
                t36 = nr*rhop*t6
                t47 = -t6
                t49 = t7-2.0_r64
                t50 = 1.0_r64/t8
                t54 = chi*t7*ve
                t55 = t10*ve*4.0_r64
                t57 = t10*t19
                t59 = t9*vk*4.0_r64
                t60 = t10*vk*4.0_r64
                t61 = t9*t21
                t64 = nps*rho*t2*2.0_r64
                t65 = nr*rhop*t3*2.0_r64
                t66 = 1.0_r64/t28
                t72 = t12*ve*8.0_r64
                t73 = t10*ve*2.1e+1_r64
                t74 = t9*ve*3.8e+1_r64
                t75 = t9*ve*4.6e+1_r64
                t76 = t9*vk*-2.0_r64
                t79 = t10*vk*-6.0_r64
                t81 = t12*vk*8.0_r64
                t83 = t9*vk*1.6e+1_r64
                t84 = t10*t33
                t86 = t9*vk*2.0e+1_r64
                t87 = t2+t5
                t88 = t18+t31
                t93 = t9*vk*-1.5e+1_r64
                t96 = 1.0_r64/t4**(3.0_r64/2.0_r64)
                t97 = 1.0_r64/t4**(5.0_r64/2.0_r64)
                t98 = 1.0_r64/t4**(7.0_r64/2.0_r64)
                t103 = t20+t37+t41
                t51 = t50**2
                t56 = t11*ve*4.0_r64
                t62 = t11*vk*4.0_r64
                t67 = t66**2
                t68 = t66**3
                t69 = -t36
                t70 = -t55
                t77 = -t59
                t78 = -t60
                t80 = t11*vk*8.0_r64
                t85 = t11*vk*1.6e+1_r64
                t89 = 1.0_r64/t34
                t91 = -t73
                t92 = -t81
                t94 = -t83
                t95 = -t86
                t99 = 1.0_r64/t49
                t102 = t3+t47
                t107 = t19+t24+t54+t76
                t114 = t32+t40+t52+t57+t61+t79
                t71 = -t56
                t90 = t89**2
                t100 = t99**2
                t101 = t99**3
                t104 = t48*t99*ve*4.0_r64
                t109 = t29+t45+t75+t94
                t110 = t48*t89*t99*(1.6e+1_r64/3.0_r64)
                t113 = t21+t38+t42+t60+t70+t77
                t115 = t48*t88*t89*t99*(8.0_r64/3.0_r64)
                t121 = t33+t39+t53+t72+t80+t84+t91+t92+t93
                t105 = t48*t100*(1.6e+1_r64/3.0_r64)
                t108 = t48*t101*ve*(1.28e+2_r64/1.5e+1_r64)
                t111 = t48*t89*t100*ve*(1.6e+1_r64/1.5e+1_r64)
                t116 = t48*t88*t90*t99*(3.2e+1_r64/1.5e+1_r64)
                t117 = t48*t88*t89*t100*(3.2e+1_r64/1.5e+1_r64)
                t118 = t30+t46+t62+t71+t74+t78+t95
                t106 = t105*ve
                t112 = -t111
                t119 = t105+t110
                t120 = t106+t115
                t122 = t88*t89*t119*(2.0_r64/5.0_r64)
                t123 = t108+t112+t116+t117+t122
                et1 = t27*t48*t51*t68*t118*(nr*t3*t17+nr*t17*t47+nps*t2*t15+nps*t5*t15+t2*t3*t4*4.0_r64-t2*t4*t6*2.0_r64+ &
                     t3*t5*t26+t4*t5*t47)*(-1.0_r64/2.0_r64)+(t27*t48*t51*t68*t114*(t2*t3*t14*2.0_r64+t2*t3*t16*2.0_r64+ &
                     t3*t5*t14*2.0_r64-t2*t6*t16*2.0_r64+t3*t5*t16+t2*t14*t47+t5*t14*t47+t5*t16*t47))/2.0_r64+ &
                     (t4*t27*t48*t51*t68*(t35+t64+t65+t69)*(t85+vk*3.0e+1_r64+chi*ve*9.0e+1_r64-chi*vk*3.0e+1_r64- &
                     t10*ve*4.2e+1_r64+t12*ve*1.6e+1_r64-t9*vk*3.0e+1_r64+t10*vk*3.0e+1_r64-t12*vk*1.6e+1_r64))/2.0_r64
                et2 = t2*t3*t4*t27*t48*t51*t68*(ve*1.5e+1_r64-chi*vk*6.0e+1_r64-t9*ve*1.35e+2_r64+t11*ve*1.36e+2_r64- &
                     t13*ve*4.8e+1_r64+t9*vk*6.0e+1_r64+t10*vk*1.0e+2_r64-t11*vk*1.0e+2_r64-t12*vk*4.8e+1_r64+ &
                     t13*vk*4.8e+1_r64)-(t4*t27*t48*t51*t68*t87*t102*t109)/2.0_r64
                et3 = -1.0_r64/t26**(5.0_r64/2.0_r64)*(t120*(t2*t3*t27*(3.0_r64/4.0_r64)-t2*t6*t27*(3.0_r64/2.0_r64)+ &
                     t3*t5*t27*(3.0_r64/2.0_r64)-t5*t6*t27*3.0_r64)+(t104-chi*t120)*(t27*t35*(3.0_r64/2.0_r64)- &
                     t27*t36*(3.0_r64/2.0_r64)+nr*nps*t22*t27*(3.0_r64/4.0_r64)+nr*rhop*t3*t27*(3.0_r64/4.0_r64)+ &
                     nps*rho*t2*t27*(3.0_r64/4.0_r64))+t2*t3*t27*(t20*t48+t9*t120-chi*t48*t99*ve*8.0_r64)*(3.0_r64/4.0_r64))
                et4 = 1.0_r64/t26**(7.0_r64/2.0_r64)*((t120-chi*t123)*(t22*t27*t35*(1.5e+1_r64/2.0_r64)-t22*t27*t36*(1.5e+1_r64/2.0_r64)+ &
                     nr*rhop*t3*t22*t27*(1.5e+1_r64/2.0_r64)+nps*rho*t2*t22*t27*(1.5e+1_r64/2.0_r64))+t123*(t2*t3*t22*t27*(1.5e+1_r64/2.0_r64)- &
                     t2*t6*t22*t27*(1.5e+1_r64/2.0_r64)+t3*t5*t22*t27*(1.5e+1_r64/2.0_r64)-t5*t6*t22*t27*(1.5e+1_r64/2.0_r64))+ &
                     t2*t3*t22*t27*(t104-chi*t120*2.0_r64+t9*t123)*(1.5e+1_r64/2.0_r64))+nz*nzs*t18*1.0_r64/t26**(3.0_r64/2.0_r64)*t27*t48*t99
                Irr = mu*((t48*t97*((t27*t48*t50*t67*t107*(t2*t3*-4.0_r64+t2*t6*2.0_r64-t3*t5*2.0_r64+t5*t6+ &
                     nz*nzs*t14+nz*nzs*t16))/4.0_r64+(t27*t48*t50*t67*t113*(t35*2.0_r64-t36*2.0_r64-nz*nzs*t4+ &
                     nr*rhop*t3*4.0_r64+nps*rho*t2*4.0_r64))/4.0_r64+t2*t3*t27*t48*t50*t67*(t19-t85-chi*vk*1.8e+1_r64- &
                     t9*ve*3.0e+1_r64+t11*ve*1.6e+1_r64+t9*vk*1.8e+1_r64+t10*vk*1.6e+1_r64)-(nz*nzs*t4*t27*t48*t50*t67*(vk+ &
                     chi*t31+chi*ve*4.0_r64))/2.0_r64))/8.0_r64+(t48*t98*(et1+et2))/1.6e+1_r64+(nr*nps*t27*t66*t96*ve)/2.0_r64)
                Irz = mu*(t48*t97*(t27*t48*t50*t67*t107*(nps*nz*t14+nz*rhop*t3-nzs*rhop*t5*2.0_r64+nps*t2*zh*2.0_r64+ &
                     nps*t5*zh)*(-1.0_r64/4.0_r64)+(t3*t27*t48*t50*t67*t113*(nz*rho+nr*zh*2.0_r64))/4.0_r64+ &
                     (nz*rho*t27*t48*t50*t67*t103*(t3-t6*2.0_r64))/4.0_r64)*(-1.0_r64/8.0_r64)-(t48*t98*((t27*t48*t51*t68*t114*zh*(nps*t2*t14+ &
                     nps*t5*t14+rhop*t2*t3*2.0_r64-rhop*t2*t6*2.0_r64+rhop*t3*t5+rhop*t5*t47))/2.0_r64-(rhop*t27*t48*t51*t68*t118*zh*(t35+ &
                     t64+t69+nr*rhop*t3))/2.0_r64+rhop*t2*t3*t27*t48*t51*t68*t121*zh-(rho*t27*t48*t51*t68*t87*t102*t109*zh)/2.0_r64))/1.6e+1_r64+ &
                     (nr*nzs*t27*t66*t96*ve)/2.0_r64)
                Izr = mu*(t48*t97*(t27*t48*t50*t67*t107*(nr*nzs*t16+nzs*rho*t2+nzs*rho*t5*2.0_r64-nr*t3*zh*2.0_r64+ &
                     nr*t6*zh)*(-1.0_r64/4.0_r64)+(t2*t27*t48*t50*t67*t113*(nzs*rhop-nps*zh*2.0_r64))/4.0_r64+ &
                     (nzs*rhop*t27*t48*t50*t67*t103*(t5+t87))/4.0_r64)*(-1.0_r64/8.0_r64)+(t48*t98*((t27*t48*t51*t68*t114*zh*(nr*t3*t16+ &
                     nr*t16*t47+rho*t2*t3*2.0_r64+rho*t3*t5*2.0_r64+rho*t2*t47+rho*t5*t47))/2.0_r64-(rho*t27*t48*t51*t68*t118*zh*(t35+ &
                     t65+t69+nps*rho*t2))/2.0_r64+rho*t2*t3*t27*t48*t51*t68*t121*zh-(rhop*t27*t48*t51*t68*t87*t102*t109*zh)/2.0_r64))/1.6e+1_r64+ &
                     (nps*nz*t27*t66*t96*ve)/2.0_r64)
                Izz = mu*(et3+et4)
              end block
            end if
            sc = rhop*skws(j)
            Gf(i,    j)   = Irr*sc; Gf(i,    p+j) = Irz*sc
            Gf(nf+i, j)   = Izr*sc; Gf(nf+i, p+j) = Izz*sc
          end do
        end do
        do ib = 1, nf
          do jp = 1, p
            B(ifr(ib),    (k-1)*p+jp)      = Gf(ib,    jp)
            B(ifr(ib),    pnsa+(k-1)*p+jp) = Gf(ib,    p+jp)
            B(nt+ifr(ib), (k-1)*p+jp)      = Gf(nf+ib, jp)
            B(nt+ifr(ib), pnsa+(k-1)*p+jp) = Gf(nf+ib, p+jp)
          end do
        end do
      end if
    end do
    A(:, 1:pnp)       = matmul(B(:, 1:pnsa),        IM)
    A(:, pnp+1:2*pnp) = matmul(B(:, pnsa+1:2*pnsa), IM)
    deallocate(tin, IM, joa, sax, sanx, saws, sawxp, saxlo, saxhi, ik, ic_, ifr, tkx, tknx, tfx, tfnx, Gc, Gf, Gfold, B)
    deallocate(cc1, cc2, cc3, cc4, cc5, asq, adq, a1q, a2q, a3q, a4q, tt_az, ww_az)
  end subroutine axissymstok_dlpn_blockmat_r64

  subroutine axissymstok_dlpn_blockmat_nmode_r64(nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, mu, A)
    integer(8),   intent(in)    :: nt, p, np, M, iside, iclosed
    complex(r64), intent(in)    :: tx(nt), tnx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: sws(p*np), tpan(np+1), mu
    complex(r64), intent(inout) :: A(3*nt, 3*np*p, M+1)
    integer(8), parameter :: Ksub = 4
    integer(8) :: qq, nso, i, j, jp, iq, ia, l, e, md, ar, b, c, qo, ne, npa, nk, nkc, jj, cols
    real(r64)  :: twopi, sumwso, sumwsc, rho, rhop, zh, rr2, chi, ws, spd, rt, zti, nrp, nzp, vk, ve, Fn, An, dFn, rn, nrt, nzt
    real(r64)  :: st1, stN, tm, denom, rlo, rhi, tgi, split
    complex(r64) :: ic, zac, zbc, cc1, cc2, cc3, cc4, cc5, Slog, Dval, Dz, Dzz, Acau, Ahy, Asup
    real(r64)    :: tglp(p), wglp(p), Dp(p,p), tglq(2*p), wglq(2*p), Dq(2*p,2*p), IP2(2*p,p)
    real(r64)    :: tc(p,np), rk(p), re2(2), Lc(p,p), IPe2(2,p), IPqc(2*p,p), wsq(2*p), wsp(p)
    complex(r64) :: Ypc(p), Ypb(p), Yq(2*p), dYq(2*p), dYp(p), nvq(2*p), wxpq(2*p), ec2(2), xc(p,np)
    real(r64),    allocatable :: tin(:)
    integer(8),   allocatable :: joa(:), ci(:), nidx(:)
    logical,      allocatable :: ikq(:), ikc(:)
    complex(r64), allocatable :: zc(:), zcn(:), znear(:), znearn(:), C1(:,:,:), C2(:,:,:), C3(:,:,:), C4(:,:,:), C5(:,:,:), Gcq(:,:), Gpc(:,:), ff(:,:)
    real(r64),    allocatable :: As(:,:), A1(:,:), A2(:,:), A3(:,:), A4(:,:)
    complex(r64), allocatable :: Ad(:,:)
    real(r64) :: pi
    qq = 2*p; nso = np*p; ic = (0.0_r64,1.0_r64); pi = acos(-1.0_r64); twopi = 2.0_r64*pi
    call gauss_r64(p,     tglp, wglp, Dp)
    call gauss_r64(2_8*p, tglq, wglq, Dq)
    call lagrange_interp_r64(p, tglp, qq, tglq, IP2)
    do j = 1, np
      do i = 1, p
        xc(i,j) = sx((j-1)*p+i)
        tc(i,j) = tpan(j) + (1.0_r64+tglp(i))/2.0_r64*(tpan(j+1)-tpan(j))
      end do
    end do
    st1 = tc(1,1); stN = tc(p,np)
    allocate(tin(np*Ksub+1+512)); ne = 0
    do j = 1, np
      do i = 0, Ksub-1
        ne = ne+1; tin(ne) = tpan(j) + (tpan(j+1)-tpan(j))*real(i,r64)/real(Ksub,r64)
      end do
    end do
    ne = ne+1; tin(ne) = tpan(np+1)
    do while (tin(2) > st1)
      split = 0.5_r64*(tin(1)+tin(2)); do i = ne, 2, -1; tin(i+1) = tin(i); end do; tin(2) = split; ne = ne+1
    end do
    do while (tin(ne-1) < stN)
      split = 0.5_r64*(tin(ne)+tin(ne-1)); tin(ne+1) = tin(ne); tin(ne) = split; ne = ne+1
    end do
    npa = ne-1
    allocate(joa(npa))
    do c = 1, npa
      tm = 0.5_r64*(tin(c)+tin(c+1)); qo = 1
      do j = 1, np
        if (tm >= tpan(j) .and. tm < tpan(j+1)) qo = j
      end do
      joa(c) = qo
    end do
    allocate(ci(nt), zc(nt), zcn(nt), nidx(nt), znear(nt), znearn(nt), ikq(nt), ikc(nt))
    allocate(C1(3*nt,3*qq,M+1), C2(3*nt,3*qq,M+1), C3(3*nt,3*qq,M+1), C4(3*nt,3*qq,M+1), C5(3*nt,3*qq,M+1))
    allocate(As(qq,nt), Ad(qq,nt), A1(qq,nt), A2(qq,nt), A3(qq,nt), A4(qq,nt), Gcq(nt,qq), Gpc(nt,p), ff(p,9))
    A = (0.0_r64,0.0_r64)
    do qo = 1, np
      cols = (qo-1)*p
      do i = 1, p
        Ypc(i) = xc(i,qo)
      end do
      sumwso = 0.0_r64
      do jp = 1, p
        sumwso = sumwso + sws(cols+jp)
      end do
      nk = 0
      do i = 1, nt
        ikq(i) = (abs(tx(i)-sxlo(qo)) + abs(tx(i)-sxhi(qo))) < 1.85_r64*sumwso
        if (ikq(i)) then; nk = nk + 1; nidx(nk) = i; znear(nk) = tx(i); znearn(nk) = tnx(i); end if
      end do
      do i = 1, nt
        if (ikq(i)) cycle
        rt = real(tx(i),r64); zti = aimag(tx(i)); nrt = real(tnx(i),r64); nzt = aimag(tnx(i))
        do jp = 1, p
          jj = cols + jp; rhop = real(sx(jj),r64); zh = zti - aimag(sx(jj)); ws = sws(jj)
          nrp = real(snx(jj),r64); nzp = aimag(snx(jj)); rho = rt
          rr2 = (rho-rhop)**2 + zh*zh; chi = 1.0_r64 + rr2/(2.0_r64*rho*rhop)
          do md = 0, M
            rn = real(md,r64); call carrier_r64(chi, md, vk, ve, Fn, An, dFn)
            block
              complex(r64) :: t2, t3, t4, t5, t7, t9, t10, t12
              complex(r64) :: t13, t15, t16, t18, t19, t21, t22, t6
              complex(r64) :: t8, t11, t14, t17, t20, t23, t27, t28
              complex(r64) :: t29, t31, t32, t34, t37, t40, t41, t42
              complex(r64) :: t43, t44, t51, t52, t53, t54, t55, t56
              complex(r64) :: t57, t58, t60, t61, t65, t66, t67, t68
              complex(r64) :: t71, t72, t73, t74, t30, t33, t35, t36
              complex(r64) :: t38, t45, t46, t47, t48, t49, t62, t63
              complex(r64) :: t64, t69, t70, t75, t76, t77, t78, t39
              complex(r64) :: et1, et2, et3, et4, t50, et5, et6, et7
              complex(r64) :: et8, et9, et10, et11, et12, et13, et14, et15
              complex(r64) :: et16, et17, et18, et19, et20, et21, et22, et23
              complex(r64) :: et24, et25, et26, et27, et28, et29, et30, et31
              complex(r64) :: et32, et33, et34, et35, et36, et37, et38, et39
              complex(r64) :: et40, et41, et42, et43, et44
              complex(r64) :: coK(9), coE(9)
              real(r64) :: nr, nz
              nr = nrt; nz = nzt
              t2 = rho*rhop
              t3 = chi+1.0_r64
              t4 = chi**2
              t5 = chi**3
              t7 = chi**5
              t9 = rn**2
              t10 = rn**3
              t12 = rho**2
              t13 = rho**3
              t15 = rhop**2
              t16 = rhop**3
              t18 = zh**2
              t19 = zh**3
              t21 = 1.0_r64/pi
              t22 = chi-1.0_r64
              t6 = t4**2
              t8 = t4**3
              t11 = t9**2
              t14 = t12**2
              t17 = t15**2
              t20 = t18**2
              t23 = 1.0_r64/t12
              t27 = t9*4.0_r64
              t28 = 1.0_r64/t3
              t29 = t4-1.0_r64
              t31 = 1.0_r64/t22**2
              t32 = 1.0_r64/t22**3
              t34 = 1.0_r64/sqrt(t2)
              t37 = nz*nzp*t19*ic
              t40 = nz*nzp*t2*zh*(2.0_r64*ic)
              t41 = chi*nz*nzp*t19*(4.0_r64*ic)
              t42 = nr*nzp*rho*t18*ic
              t43 = nrp*nz*rhop*t18*ic
              t44 = nr*nrp*t2*zh*(1.0e+1_r64*ic)
              t51 = chi*nr*nrp*t12*zh*ic
              t52 = chi*nr*nrp*t15*zh*ic
              t53 = chi*nrp*nz*rho*t18*ic
              t54 = chi*nr*nzp*rho*t18*(4.0_r64*ic)
              t55 = chi*nr*nzp*rhop*t18*ic
              t56 = chi*nrp*nz*rhop*t18*(4.0_r64*ic)
              t57 = chi*nr*nrp*t2*zh*(1.2e+1_r64*ic)
              t58 = nr*nrp*t2*t5*zh*(4.0_r64*ic)
              t60 = nr*nrp*t2*t9*zh*(4.0_r64*ic)
              t61 = nz*nzp*t9*t19*(4.0_r64*ic)
              t65 = nr*nzp*rho*t9*t18*(4.0_r64*ic)
              t66 = nrp*nz*rhop*t9*t18*(4.0_r64*ic)
              t67 = nr*nrp*t2*t4*zh*(8.0_r64*ic)
              t68 = nz*nzp*t2*t4*zh*-(2.0_r64*ic)
              t71 = chi*nr*nrp*t9*t12*zh*(4.0_r64*ic)
              t72 = chi*nr*nrp*t9*t15*zh*(4.0_r64*ic)
              t73 = chi*nrp*nz*rho*t9*t18*(4.0_r64*ic)
              t74 = chi*nr*nzp*rhop*t9*t18*(4.0_r64*ic)
              t30 = t28**2
              t33 = t27-1.0_r64
              t35 = t34**5
              t36 = t34**7
              t38 = 1.0_r64/t29
              t45 = nr*nzp*1.0_r64/t34**4*4.0_r64
              t46 = nrp*nz*1.0_r64/t34**4*4.0_r64
              t47 = nz*nzp*1.0_r64/t34**4*4.0_r64
              t48 = -t43
              t49 = -t44
              t62 = -t55
              t63 = -t56
              t64 = -t57
              t69 = -t61
              t70 = -t65
              t75 = t4*t60
              t76 = -t71
              t77 = -t72
              t78 = -t73
              t39 = t38**2
              et1 = t47-chi*nr*nrp*t14*1.0e+1_r64-chi*nr*nrp*t17*1.0e+1_r64+nr*nrp*t2*t12*2.8e+1_r64+nr*nrp*t2*t15*2.8e+1_r64 &
              +nr*nrp*t5*t14*2.0_r64+nr*nrp*t5*t17*2.0_r64-nz*nzp*t12*t18*5.0_r64-nz*nzp*t15*t18*5.0_r64-nr*nzp*t13*zh*5.0_r64 &
              +nrp*nz*t16*zh*5.0_r64-chi*nr*nrp*1.0_r64/t34**4*6.8e+1_r64+nr*nrp*t5*1.0_r64/t34**4*3.6e+1_r64 &
              -nr*nrp*t7*1.0_r64/t34**4*1.6e+1_r64-nz*nzp*t4*1.0_r64/t34**4*6.0_r64+nz*nzp*t6*1.0_r64/t34**4*2.0_r64 &
              -nz*nzp*t9*1.0_r64/t34**4*4.0_r64+chi*nr*nrp*t9*t14*4.0e+1_r64+chi*nr*nrp*t9*t17*4.0e+1_r64-chi*nz*nzp*t2*t12 &
              -chi*nz*nzp*t2*t15+chi*nz*nzp*t2*t18*1.7e+1_r64-chi*nrp*nz*t13*zh*1.0e+1_r64+chi*nr*nzp*t16*zh*1.0e+1_r64 &
              -nr*nrp*t2*t4*t12*4.0_r64
              et2 = nr*nrp*t2*t6*t12*8.0_r64-nr*nrp*t2*t4*t15*4.0_r64+nr*nrp*t2*t6*t15*8.0_r64-nr*nrp*t2*t9*t12*4.0e+1_r64 &
              -nr*nrp*t2*t9*t15*4.0e+1_r64-nr*nrp*t5*t9*t14*8.0_r64-nr*nrp*t5*t9*t17*8.0_r64+nz*nzp*t2*t5*t12+nz*nzp*t2*t5*t15 &
              -nz*nzp*t2*t5*t18-nz*nzp*t4*t12*t18*3.0_r64-nz*nzp*t4*t15*t18*3.0_r64+nz*nzp*t9*t12*t18*2.0e+1_r64 &
              +nz*nzp*t9*t15*t18*2.0e+1_r64+nrp*nz*rho*t2*zh*1.9e+1_r64-nr*nzp*rhop*t2*zh*1.9e+1_r64-nr*nzp*t4*t13*zh*3.0_r64 &
              +nrp*nz*t5*t13*zh*2.0_r64+nrp*nz*t4*t16*zh*3.0_r64-nr*nzp*t5*t16*zh*2.0_r64+nr*nzp*t9*t13*zh*2.0e+1_r64 &
              -nrp*nz*t9*t16*zh*2.0e+1_r64+chi*nr*nrp*t9*1.0_r64/t34**4*1.76e+2_r64+nr*nrp*t5*t9*1.0_r64/t34**4*4.8e+1_r64
              et3 = nr*nrp*t7*t9*1.0_r64/t34**4*-3.2e+1_r64+nz*nzp*t6*t27*1.0_r64/t34**4+chi*nr*nzp*rho*t2*zh*2.6e+1_r64 &
              -chi*nrp*nz*rhop*t2*zh*2.6e+1_r64+chi*nrp*nz*t9*t13*zh*4.0e+1_r64-chi*nr*nzp*t9*t16*zh*4.0e+1_r64 &
              -nr*nrp*t2*t4*t9*t12*1.28e+2_r64+nr*nrp*t2*t6*t9*t12*4.0e+1_r64-nr*nrp*t2*t4*t9*t15*1.28e+2_r64 &
              +nr*nrp*t2*t6*t9*t15*4.0e+1_r64-nz*nzp*t2*t5*t9*t12*4.0_r64-nz*nzp*t2*t5*t9*t15*4.0_r64 &
              +nz*nzp*t4*t9*t12*t18*1.2e+1_r64+nz*nzp*t4*t9*t15*t18*1.2e+1_r64+nz*nzp*t2*t5*t18*t27+nrp*nz*rho*t2*t4*zh &
              -nr*nzp*rho*t2*t5*zh*2.0_r64+nrp*nz*rho*t2*t6*zh*4.0_r64-nrp*nz*rho*t2*t9*zh*4.0e+1_r64-nr*nzp*rhop*t2*t4*zh &
              +nrp*nz*rhop*t2*t5*zh*2.0_r64
              et4 = nr*nzp*rhop*t2*t6*zh*-4.0_r64+nr*nzp*rhop*t2*t9*zh*4.0e+1_r64+nr*nzp*t4*t9*t13*zh*1.2e+1_r64 &
              -nrp*nz*t5*t9*t13*zh*8.0_r64-nrp*nz*t4*t9*t16*zh*1.2e+1_r64+nr*nzp*t5*t9*t16*zh*8.0_r64 &
              -chi*nz*nzp*t2*t9*t18*6.8e+1_r64+chi*nz*nzp*t2*t12*t27+chi*nz*nzp*t2*t15*t27-chi*nr*nzp*rho*t2*t9*zh*1.04e+2_r64 &
              +chi*nrp*nz*rhop*t2*t9*zh*1.04e+2_r64-nrp*nz*rho*t2*t4*t9*zh*7.6e+1_r64+nr*nzp*rho*t2*t5*t9*zh*8.0_r64 &
              +nrp*nz*rho*t2*t6*t9*zh*2.0e+1_r64+nr*nzp*rhop*t2*t4*t9*zh*7.6e+1_r64-nrp*nz*rhop*t2*t5*t9*zh*8.0_r64 &
              -nr*nzp*rhop*t2*t6*t9*zh*2.0e+1_r64
              coK(1) = rhop*t21*t36*t39*(et1+et2+et3+et4)*(-1.0_r64/8.0_r64)
              t50 = 1.0_r64/t33
              et5 = nr*nrp*1.0_r64/t34**4*5.0e+1_r64+nr*nrp*t14*1.5e+1_r64+nr*nrp*t17*1.5e+1_r64+nr*nrp*t4*t14*1.9e+1_r64 &
              -nr*nrp*t6*t14*2.0_r64+nr*nrp*t4*t17*1.9e+1_r64-nr*nrp*t6*t17*2.0_r64-nr*nrp*t9*t14*6.0e+1_r64 &
              -nr*nrp*t9*t17*6.0e+1_r64+nz*nzp*t2*t12*3.0_r64+nz*nzp*t2*t15*3.0_r64-nz*nzp*t2*t18*2.1e+1_r64 &
              +nrp*nz*t13*zh*1.5e+1_r64-nr*nzp*t16*zh*1.5e+1_r64-chi*nz*nzp*1.0_r64/t34**4*1.0e+1_r64 &
              +nr*nrp*t4*1.0_r64/t34**4*1.78e+2_r64-nr*nrp*t6*1.0_r64/t34**4*5.2e+1_r64+nr*nrp*t8*1.0_r64/t34**4*1.6e+1_r64 &
              -nr*nrp*t9*1.0_r64/t34**4*2.16e+2_r64+nr*nrp*t11*1.0_r64/t34**4*1.6e+1_r64+nz*nzp*t5*1.0_r64/t34**4*1.2e+1_r64 &
              -nz*nzp*t7*1.0_r64/t34**4*2.0_r64-chi*nr*nrp*t2*t12*1.24e+2_r64
              et6 = chi*nr*nrp*t2*t15*-1.24e+2_r64+chi*nz*nzp*t12*t18*2.9e+1_r64+chi*nz*nzp*t15*t18*2.9e+1_r64 &
              +chi*nr*nzp*t13*zh*2.9e+1_r64-chi*nrp*nz*t16*zh*2.9e+1_r64+nr*nrp*t2*t5*t12*4.0_r64-nr*nrp*t2*t7*t12*8.0_r64 &
              +nr*nrp*t2*t5*t15*4.0_r64-nr*nrp*t2*t7*t15*8.0_r64-nr*nrp*t4*t9*t14*8.0e+1_r64+nr*nrp*t4*t11*t14*1.6e+1_r64 &
              +nr*nrp*t6*t9*t14*1.2e+1_r64-nr*nrp*t4*t9*t17*8.0e+1_r64-nr*nrp*t6*t11*t14*1.6e+1_r64 &
              +nr*nrp*t4*t11*t17*1.6e+1_r64+nr*nrp*t6*t9*t17*1.2e+1_r64-nr*nrp*t6*t11*t17*1.6e+1_r64-nz*nzp*t2*t4*t12*2.0_r64 &
              -nz*nzp*t2*t6*t12-nz*nzp*t2*t4*t15*2.0_r64-nz*nzp*t2*t6*t15-nz*nzp*t2*t9*t12*1.2e+1_r64 &
              -nz*nzp*t2*t4*t18*4.4e+1_r64
              et7 = nz*nzp*t2*t6*t18-nz*nzp*t2*t9*t15*1.2e+1_r64+nz*nzp*t2*t9*t18*8.8e+1_r64-nz*nzp*t2*t11*t18*1.6e+1_r64 &
              +nz*nzp*t5*t12*t18*3.0_r64+nz*nzp*t5*t15*t18*3.0_r64-nr*nzp*rho*t2*zh*3.3e+1_r64+nrp*nz*rhop*t2*zh*3.3e+1_r64 &
              +nrp*nz*t4*t13*zh*1.9e+1_r64+nr*nzp*t5*t13*zh*3.0_r64-nrp*nz*t6*t13*zh*2.0_r64-nr*nzp*t4*t16*zh*1.9e+1_r64 &
              -nrp*nz*t5*t16*zh*3.0_r64+nr*nzp*t6*t16*zh*2.0_r64-nrp*nz*t9*t13*zh*6.0e+1_r64+nr*nzp*t9*t16*zh*6.0e+1_r64 &
              +chi*nz*nzp*t9*1.0_r64/t34**4*4.0e+1_r64-nr*nrp*t4*t9*1.0_r64/t34**4*6.88e+2_r64 &
              +nr*nrp*t4*t11*1.0_r64/t34**4*4.8e+1_r64+nr*nrp*t6*t9*1.0_r64/t34**4*1.84e+2_r64 &
              -nr*nrp*t6*t11*1.0_r64/t34**4*4.8e+1_r64-nr*nrp*t8*t9*1.0_r64/t34**4*4.8e+1_r64
              et8 = nr*nrp*t8*t11*1.0_r64/t34**4*-1.6e+1_r64-nz*nzp*t5*t9*1.0_r64/t34**4*4.8e+1_r64 &
              +nz*nzp*t7*t9*1.0_r64/t34**4*8.0_r64-chi*nrp*nz*rho*t2*zh*9.1e+1_r64+chi*nr*nzp*rhop*t2*zh*9.1e+1_r64 &
              -chi*nr*nzp*t9*t13*zh*1.2e+2_r64+chi*nr*nzp*t11*t13*zh*1.6e+1_r64+chi*nrp*nz*t9*t16*zh*1.2e+2_r64 &
              -chi*nrp*nz*t11*t16*zh*1.6e+1_r64-nr*nrp*t2*t5*t9*t12*1.6e+1_r64+nr*nrp*t2*t7*t9*t12*2.4e+1_r64 &
              -nr*nrp*t2*t5*t9*t15*1.6e+1_r64+nr*nrp*t2*t7*t11*t12*3.2e+1_r64+nr*nrp*t2*t7*t9*t15*2.4e+1_r64 &
              +nr*nrp*t2*t7*t11*t15*3.2e+1_r64+nz*nzp*t2*t4*t9*t12*8.0_r64+nz*nzp*t2*t4*t9*t15*8.0_r64 &
              +nz*nzp*t2*t4*t9*t18*1.76e+2_r64-nz*nzp*t2*t6*t9*t18*8.0_r64+nz*nzp*t2*t6*t11*t18*1.6e+1_r64
              et9 = nz*nzp*t5*t9*t12*t18*-8.0_r64-nz*nzp*t5*t11*t12*t18*1.6e+1_r64+nz*nzp*t2*t6*t12*t27 &
              -nz*nzp*t5*t9*t15*t18*8.0_r64-nz*nzp*t5*t11*t15*t18*1.6e+1_r64+nz*nzp*t2*t6*t15*t27 &
              -nr*nzp*rho*t2*t4*zh*6.5e+1_r64-nrp*nz*rho*t2*t5*zh+nr*nzp*rho*t2*t6*zh*2.0_r64-nrp*nz*rho*t2*t7*zh*4.0_r64 &
              +nr*nzp*rho*t2*t9*zh*1.36e+2_r64-nr*nzp*rho*t2*t11*zh*1.6e+1_r64+nrp*nz*rhop*t2*t4*zh*6.5e+1_r64 &
              +nr*nzp*rhop*t2*t5*zh-nrp*nz*rhop*t2*t6*zh*2.0_r64+nr*nzp*rhop*t2*t7*zh*4.0_r64-nrp*nz*rhop*t2*t9*zh*1.36e+2_r64 &
              +nrp*nz*rhop*t2*t11*zh*1.6e+1_r64-nrp*nz*t4*t9*t13*zh*8.0e+1_r64-nr*nzp*t5*t9*t13*zh*8.0_r64 &
              +nrp*nz*t4*t11*t13*zh*1.6e+1_r64
              et10 = nrp*nz*t6*t9*t13*zh*1.2e+1_r64+nr*nzp*t4*t9*t16*zh*8.0e+1_r64-nr*nzp*t5*t11*t13*zh*1.6e+1_r64 &
              +nrp*nz*t5*t9*t16*zh*8.0_r64-nrp*nz*t6*t11*t13*zh*1.6e+1_r64-nr*nzp*t4*t11*t16*zh*1.6e+1_r64 &
              -nr*nzp*t6*t9*t16*zh*1.2e+1_r64+nrp*nz*t5*t11*t16*zh*1.6e+1_r64+nr*nzp*t6*t11*t16*zh*1.6e+1_r64 &
              +chi*nr*nrp*t2*t9*t12*5.04e+2_r64-chi*nr*nrp*t2*t11*t12*3.2e+1_r64+chi*nr*nrp*t2*t9*t15*5.04e+2_r64 &
              -chi*nr*nrp*t2*t11*t15*3.2e+1_r64-chi*nz*nzp*t9*t12*t18*1.2e+2_r64+chi*nz*nzp*t11*t12*t18*1.6e+1_r64 &
              -chi*nz*nzp*t9*t15*t18*1.2e+2_r64+chi*nz*nzp*t11*t15*t18*1.6e+1_r64+chi*nrp*nz*rho*t2*t9*zh*3.72e+2_r64 &
              -chi*nrp*nz*rho*t2*t11*zh*3.2e+1_r64-chi*nr*nzp*rhop*t2*t9*zh*3.72e+2_r64
              et11 = chi*nr*nzp*rhop*t2*t11*zh*3.2e+1_r64+nr*nzp*rho*t2*t4*t9*zh*2.64e+2_r64-nr*nzp*rho*t2*t4*t11*zh*1.6e+1_r64 &
              -nr*nzp*rho*t2*t6*t9*zh*1.6e+1_r64+nrp*nz*rho*t2*t5*t11*zh*1.6e+1_r64+nrp*nz*rho*t2*t7*t9*zh*1.2e+1_r64 &
              +nr*nzp*rho*t2*t6*t11*zh*3.2e+1_r64+nrp*nz*rho*t2*t7*t11*zh*1.6e+1_r64-nrp*nz*rhop*t2*t4*t9*zh*2.64e+2_r64 &
              +nrp*nz*rhop*t2*t4*t11*zh*1.6e+1_r64+nrp*nz*rhop*t2*t6*t9*zh*1.6e+1_r64-nr*nzp*rhop*t2*t5*t11*zh*1.6e+1_r64 &
              -nr*nzp*rhop*t2*t7*t9*zh*1.2e+1_r64-nrp*nz*rhop*t2*t6*t11*zh*3.2e+1_r64-nr*nzp*rhop*t2*t7*t11*zh*1.6e+1_r64
              coE(1) = (t21*t30*t32*t35*t50*(et5+et6+et7+et8+et9+et10+et11))/(rho*8.0_r64)
              et12 = nr*nrp*t16*-(9.0_r64*ic)+chi*nr*nrp*t13*ic-nr*nrp*rho*t2*(1.9e+1_r64*ic)+nz*nzp*rho*t18*ic &
              -nz*nzp*rhop*t2*(3.0_r64*ic)+nr*nrp*t4*t16*(8.0_r64*ic)-nrp*nz*t2*zh*(1.0e+1_r64*ic)+nr*nzp*t12*zh*ic &
              +nr*nzp*t15*zh*(9.0_r64*ic)+chi*nr*nrp*rhop*t2*(2.3e+1_r64*ic)-chi*nz*nzp*rhop*t18*ic &
              -chi*nr*nrp*t9*t13*(4.0_r64*ic)-chi*nr*nzp*t2*zh*(2.0_r64*ic)+chi*nrp*nz*t12*zh*ic+chi*nrp*nz*t15*zh*ic &
              +nr*nrp*rho*t2*t4*(1.6e+1_r64*ic)+nr*nrp*rho*t2*t9*(4.0_r64*ic)-nr*nrp*rhop*t2*t5*(2.0e+1_r64*ic) &
              -nz*nzp*rho*t9*t18*(4.0_r64*ic)+nz*nzp*rhop*t2*t4*(3.0_r64*ic)+nr*nrp*t4*t9*t16*(4.0_r64*ic)
              et13 = nrp*nz*t2*t4*zh*(8.0_r64*ic)+nrp*nz*t2*t9*zh*(4.0_r64*ic)-nr*nzp*t4*t15*zh*(8.0_r64*ic) &
              -nr*nzp*t9*t12*zh*(4.0_r64*ic)+chi*nr*nzp*t2*t9*zh*(8.0_r64*ic)-chi*nrp*nz*t9*t12*zh*(4.0_r64*ic) &
              -chi*nrp*nz*t9*t15*zh*(4.0_r64*ic)+nr*nrp*rho*t2*t4*t9*(8.0_r64*ic)-nr*nrp*rhop*t2*t5*t9*(4.0_r64*ic) &
              +nrp*nz*t2*t4*t9*zh*(4.0_r64*ic)-nr*nzp*t4*t9*t15*zh*(4.0_r64*ic)-chi*nr*nrp*rhop*t2*t9*(8.0_r64*ic) &
              +chi*nz*nzp*rhop*t9*t18*(4.0_r64*ic)
              coK(2) = rn*t21*t35*t38*(et12+et13)*(-1.0_r64/4.0_r64)
              et14 = rn*nr*nrp*t13*(5.0_r64*ic)-nr*nrp*t10*t13*(2.0e+1_r64*ic)-chi*rn*nr*nrp*t16*(8.0_r64*ic) &
              +chi*nr*nrp*t10*t16*(3.2e+1_r64*ic)+rn*nr*nrp*rhop*t2*(1.1e+1_r64*ic)+rn*nz*nzp*rho*t2*ic &
              -rn*nz*nzp*rhop*t18*(4.0_r64*ic)-rn*nr*nrp*t4*t13*ic+rn*nr*nrp*t5*t16*(4.0_r64*ic) &
              -nr*nrp*rhop*t2*t10*(3.2e+1_r64*ic)-nz*nzp*rho*t2*t10*(4.0_r64*ic)+nz*nzp*rhop*t10*t18*(1.6e+1_r64*ic) &
              +nr*nrp*t4*t10*t13*(4.0_r64*ic)-nr*nrp*t5*t10*t16*(1.6e+1_r64*ic)-rn*nr*nzp*t2*zh*(8.0_r64*ic) &
              +rn*nrp*nz*t12*zh*(5.0_r64*ic)+rn*nrp*nz*t15*zh*(4.0_r64*ic)+nr*nzp*t2*t10*zh*(3.2e+1_r64*ic) &
              -nrp*nz*t10*t12*zh*(2.0e+1_r64*ic)-nrp*nz*t10*t15*zh*(1.6e+1_r64*ic)
              et15 = chi*nrp*nz*t2*t10*zh*(4.8e+1_r64*ic)-chi*nr*nzp*t10*t12*zh*(1.6e+1_r64*ic) &
              -chi*nr*nzp*t10*t15*zh*(3.2e+1_r64*ic)+rn*nr*nrp*rho*t2*t5*(8.0_r64*ic)+rn*nr*nrp*rhop*t2*t4*(5.0_r64*ic) &
              -rn*nr*nrp*rhop*t2*t6*(4.0_r64*ic)-rn*nz*nzp*rho*t2*t4*ic+rn*nz*nzp*rhop*t2*t5*ic &
              -nr*nrp*rho*t2*t5*t10*(3.2e+1_r64*ic)-nr*nrp*rhop*t2*t4*t10*(4.4e+1_r64*ic)+nr*nrp*rhop*t2*t6*t10*(2.8e+1_r64*ic) &
              +nz*nzp*rho*t2*t4*t10*(4.0_r64*ic)-nz*nzp*rhop*t2*t5*t10*(4.0_r64*ic)+rn*nrp*nz*t2*t5*zh*(4.0_r64*ic) &
              -rn*nrp*nz*t4*t12*zh*ic-rn*nr*nzp*t5*t15*zh*(4.0_r64*ic)-nrp*nz*t2*t5*t10*zh*(1.6e+1_r64*ic) &
              +nrp*nz*t4*t10*t12*zh*(4.0_r64*ic)
              et16 = nr*nzp*t5*t10*t15*zh*(1.6e+1_r64*ic)-chi*rn*nr*nrp*rho*t2*(2.0e+1_r64*ic) &
              +chi*rn*nz*nzp*rho*t18*(4.0_r64*ic)-chi*rn*nz*nzp*rhop*t2*ic+chi*nr*nrp*rho*t2*t10*(8.0e+1_r64*ic) &
              -chi*nz*nzp*rho*t10*t18*(1.6e+1_r64*ic)+chi*nz*nzp*rhop*t2*t10*(4.0_r64*ic)-chi*rn*nrp*nz*t2*zh*(1.2e+1_r64*ic) &
              +chi*rn*nr*nzp*t12*zh*(4.0_r64*ic)+chi*rn*nr*nzp*t15*zh*(8.0_r64*ic)
              coE(2) = t21*t28*t31*t35*t50*(et14+et15+et16)*(-1.0_r64/4.0_r64)
              et17 = -t46-nz*nzp*rhop*t19*5.0_r64-nr*nzp*t2*t18*1.0e+1_r64+nrp*nz*t12*t18*5.0_r64+nrp*nz*t15*t18*5.0_r64 &
              +nr*nrp*t13*zh*5.0_r64+nrp*nz*t4*1.0_r64/t34**4*6.0_r64-nrp*nz*t6*1.0_r64/t34**4*2.0_r64 &
              +nrp*nz*t27*1.0_r64/t34**4+chi*nz*nzp*rho*t19*8.0_r64+chi*nrp*nz*t2*t12+chi*nrp*nz*t2*t15 &
              -chi*nrp*nz*t2*t18*1.7e+1_r64+chi*nr*nzp*t12*t18*8.0_r64+chi*nr*nzp*t15*t18*1.0e+1_r64 &
              -chi*nr*nrp*t16*zh*1.0e+1_r64-nz*nzp*rhop*t4*t19*3.0_r64+nz*nzp*rhop*t9*t19*2.0e+1_r64-nrp*nz*t2*t5*t12 &
              -nrp*nz*t2*t5*t15-nr*nzp*t2*t4*t18*6.0_r64+nrp*nz*t2*t5*t18+nr*nzp*t2*t9*t18*4.0e+1_r64+nrp*nz*t4*t12*t18*3.0_r64 &
              +nrp*nz*t4*t15*t18*3.0_r64-nr*nzp*t5*t15*t18*2.0_r64
              et18 = nrp*nz*t9*t12*t18*-2.0e+1_r64-nrp*nz*t9*t15*t18*2.0e+1_r64+nr*nrp*rhop*t2*zh*1.9e+1_r64 &
              +nz*nzp*rho*t2*zh*2.0_r64+nr*nrp*t4*t13*zh*3.0_r64+nr*nrp*t5*t16*zh*2.0_r64-nr*nrp*t9*t13*zh*2.0e+1_r64 &
              -nrp*nz*t6*t9*1.0_r64/t34**4*4.0_r64-chi*nr*nrp*rho*t2*zh*2.6e+1_r64-chi*nz*nzp*rhop*t2*zh*2.0_r64 &
              +chi*nr*nrp*t9*t16*zh*4.0e+1_r64+nz*nzp*rhop*t4*t9*t19*1.2e+1_r64+nr*nzp*t2*t4*t9*t18*2.4e+1_r64 &
              -nrp*nz*t2*t5*t9*t18*4.0_r64-nrp*nz*t4*t9*t12*t18*1.2e+1_r64+nrp*nz*t2*t5*t12*t27-nrp*nz*t4*t9*t15*t18*1.2e+1_r64 &
              +nr*nzp*t5*t9*t15*t18*8.0_r64+nrp*nz*t2*t5*t15*t27+nr*nrp*rho*t2*t5*zh*2.0_r64+nr*nrp*rhop*t2*t4*zh &
              +nr*nrp*rhop*t2*t6*zh*4.0_r64-nr*nrp*rhop*t2*t9*zh*4.0e+1_r64
              et19 = nz*nzp*rho*t2*t4*zh*-2.0_r64-nz*nzp*rho*t2*t9*zh*8.0_r64+nz*nzp*rhop*t2*t5*zh*2.0_r64 &
              -nr*nrp*t4*t9*t13*zh*1.2e+1_r64-nr*nrp*t5*t9*t16*zh*8.0_r64-chi*nz*nzp*rho*t9*t19*3.2e+1_r64 &
              -chi*nrp*nz*t2*t9*t12*4.0_r64-chi*nrp*nz*t2*t9*t15*4.0_r64+chi*nrp*nz*t2*t9*t18*6.8e+1_r64 &
              -chi*nr*nzp*t9*t12*t18*3.2e+1_r64-chi*nr*nzp*t9*t15*t18*4.0e+1_r64+chi*nr*nrp*rho*t2*t9*zh*1.04e+2_r64 &
              +chi*nz*nzp*rhop*t2*t9*zh*8.0_r64-nr*nrp*rho*t2*t5*t9*zh*8.0_r64-nr*nrp*rhop*t2*t4*t9*zh*7.6e+1_r64 &
              +nr*nrp*rhop*t2*t6*t9*zh*2.0e+1_r64+nz*nzp*rho*t2*t4*t9*zh*8.0_r64-nz*nzp*rhop*t2*t5*t9*zh*8.0_r64
              coK(3) = (rhop*t21*t36*t39*(et17+et18+et19))/8.0_r64
              et20 = -t45+nz*nzp*rho*t19*9.0_r64+nrp*nz*t2*t12*3.0_r64+nrp*nz*t2*t15*3.0_r64-nrp*nz*t2*t18*2.1e+1_r64 &
              +nr*nzp*t12*t18*9.0_r64+nr*nzp*t15*t18*1.5e+1_r64-nr*nrp*t16*zh*1.5e+1_r64-chi*nrp*nz*1.0_r64/t34**4*1.0e+1_r64 &
              +nr*nzp*t4*1.0_r64/t34**4*8.0_r64+nrp*nz*t5*1.0_r64/t34**4*1.2e+1_r64-nr*nzp*t6*1.0_r64/t34**4*4.0_r64 &
              -nrp*nz*t7*1.0_r64/t34**4*2.0_r64-chi*nz*nzp*rhop*t19*2.9e+1_r64-chi*nr*nzp*t2*t18*5.8e+1_r64 &
              +chi*nrp*nz*t12*t18*2.9e+1_r64+chi*nrp*nz*t15*t18*2.9e+1_r64+chi*nr*nrp*t13*zh*2.9e+1_r64 &
              +nz*nzp*rho*t4*t19*2.3e+1_r64-nz*nzp*rho*t9*t19*4.0_r64-nz*nzp*rhop*t5*t19*3.0_r64-nrp*nz*t2*t4*t12*2.0_r64 &
              -nrp*nz*t2*t6*t12-nrp*nz*t2*t4*t15*2.0_r64
              et21 = -nrp*nz*t2*t6*t15-nrp*nz*t2*t4*t18*4.4e+1_r64-nr*nzp*t2*t5*t18*6.0_r64+nrp*nz*t2*t6*t18 &
              +nr*nzp*t4*t12*t18*2.3e+1_r64+nrp*nz*t5*t12*t18*3.0_r64+nr*nzp*t4*t15*t18*1.9e+1_r64+nrp*nz*t5*t15*t18*3.0_r64 &
              -nr*nzp*t6*t15*t18*2.0_r64-nr*nzp*t9*t12*t18*4.0_r64+nrp*nz*t2*t18*t27-nr*nrp*rho*t2*zh*3.3e+1_r64 &
              -nz*nzp*rhop*t2*zh*6.0_r64+nr*nrp*t5*t13*zh*3.0_r64-nr*nrp*t4*t16*zh*1.9e+1_r64+nr*nrp*t6*t16*zh*2.0_r64 &
              +chi*nr*nrp*rhop*t2*zh*9.1e+1_r64+chi*nz*nzp*rho*t2*zh*8.0_r64-chi*nr*nrp*t9*t13*zh*4.0_r64+nz*nzp*rho*t4*t19*t27 &
              -nz*nzp*rhop*t5*t9*t19*4.0_r64-nr*nzp*t2*t5*t9*t18*8.0_r64-nrp*nz*t2*t6*t9*t18*4.0_r64 &
              -nr*nzp*t4*t9*t15*t18*4.0_r64
              et22 = nr*nzp*t4*t12*t18*t27+nrp*nz*t5*t12*t18*t27+nrp*nz*t5*t15*t18*t27+nr*nzp*t6*t15*t18*t27 &
              -nr*nrp*rho*t2*t4*zh*6.5e+1_r64+nr*nrp*rho*t2*t6*zh*2.0_r64+nr*nrp*rho*t2*t27*zh+nr*nrp*rhop*t2*t5*zh &
              +nr*nrp*rhop*t2*t7*zh*4.0_r64-nz*nzp*rho*t2*t5*zh*8.0_r64+nz*nzp*rhop*t2*t4*zh*4.0_r64 &
              +nz*nzp*rhop*t2*t6*zh*2.0_r64-nr*nrp*t6*t9*t16*zh*4.0_r64+nr*nrp*t5*t13*t27*zh+nr*nrp*t4*t16*t27*zh &
              +chi*nz*nzp*rhop*t19*t27+chi*nr*nzp*t2*t9*t18*8.0_r64-chi*nrp*nz*t9*t12*t18*4.0_r64-chi*nrp*nz*t9*t15*t18*4.0_r64 &
              -chi*nr*nrp*rhop*t2*t9*zh*8.0_r64-nr*nrp*rho*t2*t6*t9*zh*8.0_r64+nr*nrp*rho*t2*t4*t27*zh+nr*nrp*rhop*t2*t5*t27*zh &
              +nr*nrp*rhop*t2*t7*t27*zh
              coE(3) = rhop*t21*t30*t32*t36*(et20+et21+et22)*(-1.0_r64/8.0_r64)
              et23 = nr*nrp*t13*(9.0_r64*ic)-chi*nr*nrp*t16*ic+nr*nrp*rhop*t2*(1.9e+1_r64*ic)+nz*nzp*rho*t2*(3.0_r64*ic) &
              -nz*nzp*rhop*t18*ic-nr*nrp*t4*t13*(8.0_r64*ic)-nr*nzp*t2*zh*(1.0e+1_r64*ic)+nrp*nz*t12*zh*(9.0_r64*ic) &
              +nrp*nz*t15*zh*ic-chi*nr*nrp*rho*t2*(2.3e+1_r64*ic)+chi*nz*nzp*rho*t18*ic+chi*nr*nrp*t9*t16*(4.0_r64*ic) &
              -chi*nrp*nz*t2*zh*(2.0_r64*ic)+chi*nr*nzp*t12*zh*ic+chi*nr*nzp*t15*zh*ic+nr*nrp*rho*t2*t5*(2.0e+1_r64*ic) &
              -nr*nrp*rhop*t2*t4*(1.6e+1_r64*ic)-nr*nrp*rhop*t2*t9*(4.0_r64*ic)-nz*nzp*rho*t2*t4*(3.0_r64*ic) &
              +nz*nzp*rhop*t9*t18*(4.0_r64*ic)-nr*nrp*t4*t9*t13*(4.0_r64*ic)
              et24 = nr*nzp*t2*t4*zh*(8.0_r64*ic)+nr*nzp*t2*t9*zh*(4.0_r64*ic)-nrp*nz*t4*t12*zh*(8.0_r64*ic) &
              -nrp*nz*t9*t15*zh*(4.0_r64*ic)+chi*nrp*nz*t2*t9*zh*(8.0_r64*ic)-chi*nr*nzp*t9*t12*zh*(4.0_r64*ic) &
              -chi*nr*nzp*t9*t15*zh*(4.0_r64*ic)+nr*nrp*rho*t2*t5*t9*(4.0_r64*ic)-nr*nrp*rhop*t2*t4*t9*(8.0_r64*ic) &
              +nr*nzp*t2*t4*t9*zh*(4.0_r64*ic)-nrp*nz*t4*t9*t12*zh*(4.0_r64*ic)+chi*nr*nrp*rho*t2*t9*(8.0_r64*ic) &
              -chi*nz*nzp*rho*t9*t18*(4.0_r64*ic)
              coK(4) = rn*t15*t21*t36*t38*(et23+et24)*(-1.0_r64/4.0_r64)
              et25 = rn*nr*nrp*t16*(5.0_r64*ic)-nr*nrp*t10*t16*(2.0e+1_r64*ic)-chi*rn*nr*nrp*t13*(8.0_r64*ic) &
              +chi*nr*nrp*t10*t13*(3.2e+1_r64*ic)+rn*nr*nrp*rho*t2*(1.1e+1_r64*ic)-rn*nz*nzp*rho*t18*(4.0_r64*ic) &
              +rn*nz*nzp*rhop*t2*ic+rn*nr*nrp*t5*t13*(4.0_r64*ic)-rn*nr*nrp*t4*t16*ic-nr*nrp*rho*t2*t10*(3.2e+1_r64*ic) &
              +nz*nzp*rho*t10*t18*(1.6e+1_r64*ic)-nz*nzp*rhop*t2*t10*(4.0_r64*ic)-nr*nrp*t5*t10*t13*(1.6e+1_r64*ic) &
              +nr*nrp*t4*t10*t16*(4.0_r64*ic)+rn*nrp*nz*t2*zh*(8.0_r64*ic)-rn*nr*nzp*t12*zh*(4.0_r64*ic) &
              -rn*nr*nzp*t15*zh*(5.0_r64*ic)-nrp*nz*t2*t10*zh*(3.2e+1_r64*ic)+nr*nzp*t10*t12*zh*(1.6e+1_r64*ic) &
              +nr*nzp*t10*t15*zh*(2.0e+1_r64*ic)
              et26 = chi*nr*nzp*t2*t10*zh*-(4.8e+1_r64*ic)+chi*nrp*nz*t10*t12*zh*(3.2e+1_r64*ic) &
              +chi*nrp*nz*t10*t15*zh*(1.6e+1_r64*ic)+rn*nr*nrp*rho*t2*t4*(5.0_r64*ic)-rn*nr*nrp*rho*t2*t6*(4.0_r64*ic) &
              +rn*nr*nrp*rhop*t2*t5*(8.0_r64*ic)+rn*nz*nzp*rho*t2*t5*ic-rn*nz*nzp*rhop*t2*t4*ic &
              -nr*nrp*rho*t2*t4*t10*(4.4e+1_r64*ic)+nr*nrp*rho*t2*t6*t10*(2.8e+1_r64*ic)-nr*nrp*rhop*t2*t5*t10*(3.2e+1_r64*ic) &
              -nz*nzp*rho*t2*t5*t10*(4.0_r64*ic)+nz*nzp*rhop*t2*t4*t10*(4.0_r64*ic)-rn*nr*nzp*t2*t5*zh*(4.0_r64*ic) &
              +rn*nrp*nz*t5*t12*zh*(4.0_r64*ic)+rn*nr*nzp*t4*t15*zh*ic+nr*nzp*t2*t5*t10*zh*(1.6e+1_r64*ic) &
              -nrp*nz*t5*t10*t12*zh*(1.6e+1_r64*ic)
              et27 = nr*nzp*t4*t10*t15*zh*-(4.0_r64*ic)-chi*rn*nr*nrp*rhop*t2*(2.0e+1_r64*ic)-chi*rn*nz*nzp*rho*t2*ic &
              +chi*rn*nz*nzp*rhop*t18*(4.0_r64*ic)+chi*nr*nrp*rhop*t2*t10*(8.0e+1_r64*ic)+chi*nz*nzp*rho*t2*t10*(4.0_r64*ic) &
              -chi*nz*nzp*rhop*t10*t18*(1.6e+1_r64*ic)+chi*rn*nr*nzp*t2*zh*(1.2e+1_r64*ic)-chi*rn*nrp*nz*t12*zh*(8.0_r64*ic) &
              -chi*rn*nrp*nz*t15*zh*(4.0_r64*ic)
              coE(4) = (t21*t23*t28*t31*t34**3*t50*(et25+et26+et27))/4.0_r64
              et28 = nr*nrp*t12*5.0_r64+nr*nrp*t15*5.0_r64+nz*nzp*t2*2.0_r64-chi*nr*nrp*t2*1.8e+1_r64+chi*nz*nzp*t18 &
              +nr*nrp*t2*t5*1.6e+1_r64-nr*nrp*t4*t12*4.0_r64-nr*nrp*t4*t15*4.0_r64+nr*nrp*t9*t12*1.6e+1_r64 &
              +nr*nrp*t9*t15*1.6e+1_r64-nz*nzp*t2*t4*2.0_r64+nz*nzp*t2*t27+nrp*nz*rho*zh*5.0_r64-nr*nzp*rhop*zh*5.0_r64 &
              -chi*nr*nrp*t2*t9*2.4e+1_r64-chi*nz*nzp*t9*t18*4.0_r64+chi*nr*nzp*rho*zh-chi*nrp*nz*rhop*zh &
              +nr*nrp*t2*t5*t9*3.2e+1_r64-nr*nrp*t4*t9*t12*2.0e+1_r64-nr*nrp*t4*t9*t15*2.0e+1_r64-nz*nzp*t2*t4*t9*4.0_r64 &
              -nrp*nz*rho*t4*zh*4.0_r64+nrp*nz*rho*t9*zh*1.6e+1_r64+nr*nzp*rhop*t4*zh*4.0_r64-nr*nzp*rhop*t9*zh*1.6e+1_r64 &
              -chi*nr*nzp*rho*t9*zh*4.0_r64
              et29 = chi*nrp*nz*rhop*t27*zh-nrp*nz*rho*t4*t9*zh*2.0e+1_r64+nr*nzp*rhop*t4*t9*zh*2.0e+1_r64
              coK(5) = (t21*t34**3*t38*(et28+et29))/(rho*8.0_r64)
              et30 = nr*nrp*t2*6.0_r64+nz*nzp*t18*3.0_r64+chi*nr*nrp*t12*8.0_r64+chi*nr*nrp*t15*8.0_r64+chi*nz*nzp*t2*2.0_r64 &
              -nr*nrp*t2*t4*3.0e+1_r64+nr*nrp*t2*t6*1.6e+1_r64-nr*nrp*t2*t9*1.6e+1_r64+nr*nrp*t2*t11*1.6e+1_r64 &
              -nr*nrp*t5*t12*4.0_r64-nr*nrp*t5*t15*4.0_r64-nz*nzp*t2*t5*2.0_r64+nz*nzp*t4*t18-nz*nzp*t9*t18*8.0_r64 &
              -nz*nzp*t11*t18*1.6e+1_r64+nr*nzp*rho*zh*3.0_r64-nrp*nz*rhop*zh*3.0_r64-chi*nr*nrp*t9*t12*2.8e+1_r64 &
              -chi*nr*nrp*t11*t12*1.6e+1_r64-chi*nr*nrp*t9*t15*2.8e+1_r64-chi*nr*nrp*t11*t15*1.6e+1_r64 &
              -chi*nz*nzp*t2*t9*8.0_r64+chi*nrp*nz*rho*zh*8.0_r64-chi*nr*nzp*rhop*zh*8.0_r64+nr*nrp*t2*t4*t9*9.6e+1_r64 &
              -nr*nrp*t2*t6*t9*4.8e+1_r64
              et31 = nr*nrp*t2*t6*t11*-1.6e+1_r64+nr*nrp*t5*t9*t12*1.2e+1_r64+nr*nrp*t5*t11*t12*1.6e+1_r64 &
              +nr*nrp*t5*t9*t15*1.2e+1_r64+nr*nrp*t5*t11*t15*1.6e+1_r64+nz*nzp*t2*t5*t9*8.0_r64-nz*nzp*t4*t9*t18*8.0_r64 &
              +nz*nzp*t4*t11*t18*1.6e+1_r64+nr*nzp*rho*t4*zh-nrp*nz*rho*t5*zh*4.0_r64-nr*nzp*rho*t9*zh*8.0_r64 &
              -nr*nzp*rho*t11*zh*1.6e+1_r64-nrp*nz*rhop*t4*zh+nr*nzp*rhop*t5*zh*4.0_r64+nrp*nz*rhop*t9*zh*8.0_r64 &
              +nrp*nz*rhop*t11*zh*1.6e+1_r64-chi*nrp*nz*rho*t9*zh*2.8e+1_r64-chi*nrp*nz*rho*t11*zh*1.6e+1_r64 &
              +chi*nr*nzp*rhop*t9*zh*2.8e+1_r64+chi*nr*nzp*rhop*t11*zh*1.6e+1_r64-nr*nzp*rho*t4*t9*zh*8.0_r64 &
              +nrp*nz*rho*t5*t9*zh*1.2e+1_r64+nr*nzp*rho*t4*t11*zh*1.6e+1_r64
              et32 = nrp*nz*rho*t5*t11*zh*1.6e+1_r64+nrp*nz*rhop*t4*t9*zh*8.0_r64-nr*nzp*rhop*t5*t9*zh*1.2e+1_r64 &
              -nrp*nz*rhop*t4*t11*zh*1.6e+1_r64-nr*nzp*rhop*t5*t11*zh*1.6e+1_r64
              coE(5) = (rhop*t21*t28*t31*t35*t50*(et30+et31+et32))/8.0_r64
              coK(6) = rn*t21*t23*t34**3*t38*(t37+t42+t48+t49+t51+t52+t53+t60+t62+t66+t67+t69+t70+t74+t75+t76+t77+t78 &
              +nrp*nz*rho*t2*(3.0_r64*ic)-nrp*nz*rho*t2*t4*(3.0_r64*ic))*(-1.0_r64/4.0_r64)
              coE(6) = (rn*rhop*t21*t28*t31*t35*(t40+t41+t54+t58+t63+t64+t68+nrp*nz*rho*t18*(4.0_r64*ic)-nrp*nz*rhop*t2*ic &
              -nr*nzp*rhop*t18*(5.0_r64*ic)+nr*nrp*t12*zh*(4.0_r64*ic)+nr*nrp*t15*zh*(5.0_r64*ic)+chi*nrp*nz*rho*t2*ic &
              -nrp*nz*rho*t2*t5*ic+nrp*nz*rhop*t2*t4*ic+nr*nzp*rhop*t4*t18*ic-nr*nrp*t4*t15*zh*ic))/(rho*4.0_r64)
              et33 = t45-nz*nzp*rho*t19*5.0_r64+nrp*nz*t2*t18*1.0e+1_r64-nr*nzp*t12*t18*5.0_r64-nr*nzp*t15*t18*5.0_r64 &
              +nr*nrp*t16*zh*5.0_r64-nr*nzp*t4*1.0_r64/t34**4*6.0_r64+nr*nzp*t6*1.0_r64/t34**4*2.0_r64 &
              -nr*nzp*t9*1.0_r64/t34**4*4.0_r64+chi*nz*nzp*rhop*t19*8.0_r64-chi*nr*nzp*t2*t12-chi*nr*nzp*t2*t15 &
              +chi*nr*nzp*t2*t18*1.7e+1_r64-chi*nrp*nz*t12*t18*1.0e+1_r64-chi*nrp*nz*t15*t18*8.0_r64 &
              -chi*nr*nrp*t13*zh*1.0e+1_r64-nz*nzp*rho*t4*t19*3.0_r64+nz*nzp*rho*t9*t19*2.0e+1_r64+nr*nzp*t2*t5*t12 &
              +nr*nzp*t2*t5*t15+nrp*nz*t2*t4*t18*6.0_r64-nr*nzp*t2*t5*t18-nrp*nz*t2*t9*t18*4.0e+1_r64-nr*nzp*t4*t12*t18*3.0_r64 &
              +nrp*nz*t5*t12*t18*2.0_r64-nr*nzp*t4*t15*t18*3.0_r64
              et34 = nr*nzp*t9*t12*t18*2.0e+1_r64+nr*nzp*t9*t15*t18*2.0e+1_r64+nr*nrp*rho*t2*zh*1.9e+1_r64 &
              +nz*nzp*rhop*t2*zh*2.0_r64+nr*nrp*t5*t13*zh*2.0_r64+nr*nrp*t4*t16*zh*3.0_r64-nr*nrp*t9*t16*zh*2.0e+1_r64 &
              +nr*nzp*t6*t27*1.0_r64/t34**4-chi*nr*nrp*rhop*t2*zh*2.6e+1_r64-chi*nz*nzp*rho*t2*zh*2.0_r64 &
              +chi*nr*nrp*t9*t13*zh*4.0e+1_r64+nz*nzp*rho*t4*t9*t19*1.2e+1_r64-nr*nzp*t2*t5*t9*t12*4.0_r64 &
              -nr*nzp*t2*t5*t9*t15*4.0_r64-nrp*nz*t2*t4*t9*t18*2.4e+1_r64+nr*nzp*t4*t9*t12*t18*1.2e+1_r64 &
              -nrp*nz*t5*t9*t12*t18*8.0_r64+nr*nzp*t4*t9*t15*t18*1.2e+1_r64+nr*nzp*t2*t5*t18*t27+nr*nrp*rho*t2*t4*zh &
              +nr*nrp*rho*t2*t6*zh*4.0_r64-nr*nrp*rho*t2*t9*zh*4.0e+1_r64+nr*nrp*rhop*t2*t5*zh*2.0_r64
              et35 = nz*nzp*rho*t2*t5*zh*2.0_r64-nz*nzp*rhop*t2*t4*zh*2.0_r64-nz*nzp*rhop*t2*t9*zh*8.0_r64 &
              -nr*nrp*t5*t9*t13*zh*8.0_r64-nr*nrp*t4*t9*t16*zh*1.2e+1_r64-chi*nz*nzp*rhop*t9*t19*3.2e+1_r64 &
              -chi*nr*nzp*t2*t9*t18*6.8e+1_r64+chi*nrp*nz*t9*t12*t18*4.0e+1_r64+chi*nr*nzp*t2*t12*t27 &
              +chi*nrp*nz*t9*t15*t18*3.2e+1_r64+chi*nr*nzp*t2*t15*t27+chi*nr*nrp*rhop*t2*t9*zh*1.04e+2_r64 &
              +chi*nz*nzp*rho*t2*t9*zh*8.0_r64-nr*nrp*rho*t2*t4*t9*zh*7.6e+1_r64+nr*nrp*rho*t2*t6*t9*zh*2.0e+1_r64 &
              -nr*nrp*rhop*t2*t5*t9*zh*8.0_r64-nz*nzp*rho*t2*t5*t9*zh*8.0_r64+nz*nzp*rhop*t2*t4*t9*zh*8.0_r64
              coK(7) = rhop*t21*t36*t39*(et33+et34+et35)*(-1.0_r64/8.0_r64)
              et36 = t46+t6*t46+nz*nzp*rhop*t19*9.0_r64-nr*nzp*t2*t12*3.0_r64-nr*nzp*t2*t15*3.0_r64+nr*nzp*t2*t18*2.1e+1_r64 &
              -nrp*nz*t12*t18*1.5e+1_r64-nrp*nz*t15*t18*9.0_r64-nr*nrp*t13*zh*1.5e+1_r64+chi*nr*nzp*1.0_r64/t34**4*1.0e+1_r64 &
              -nrp*nz*t4*1.0_r64/t34**4*8.0_r64-nr*nzp*t5*1.0_r64/t34**4*1.2e+1_r64+nr*nzp*t7*1.0_r64/t34**4*2.0_r64 &
              -chi*nz*nzp*rho*t19*2.9e+1_r64+chi*nrp*nz*t2*t18*5.8e+1_r64-chi*nr*nzp*t12*t18*2.9e+1_r64 &
              -chi*nr*nzp*t15*t18*2.9e+1_r64+chi*nr*nrp*t16*zh*2.9e+1_r64-nz*nzp*rho*t5*t19*3.0_r64 &
              +nz*nzp*rhop*t4*t19*2.3e+1_r64-nz*nzp*rhop*t9*t19*4.0_r64+nr*nzp*t2*t4*t12*2.0_r64+nr*nzp*t2*t6*t12 &
              +nr*nzp*t2*t4*t15*2.0_r64+nr*nzp*t2*t6*t15+nr*nzp*t2*t4*t18*4.4e+1_r64
              et37 = nrp*nz*t2*t5*t18*6.0_r64-nr*nzp*t2*t6*t18-nr*nzp*t2*t9*t18*4.0_r64-nrp*nz*t4*t12*t18*1.9e+1_r64 &
              -nr*nzp*t5*t12*t18*3.0_r64+nrp*nz*t6*t12*t18*2.0_r64-nrp*nz*t4*t15*t18*2.3e+1_r64-nr*nzp*t5*t15*t18*3.0_r64 &
              +nrp*nz*t15*t18*t27-nr*nrp*rhop*t2*zh*3.3e+1_r64-nz*nzp*rho*t2*zh*6.0_r64-nr*nrp*t4*t13*zh*1.9e+1_r64 &
              +nr*nrp*t6*t13*zh*2.0_r64+nr*nrp*t5*t16*zh*3.0_r64+chi*nr*nrp*rho*t2*zh*9.1e+1_r64+chi*nz*nzp*rhop*t2*zh*8.0_r64 &
              -chi*nr*nrp*t9*t16*zh*4.0_r64-nz*nzp*rho*t5*t9*t19*4.0_r64+nz*nzp*rhop*t4*t19*t27+nrp*nz*t2*t5*t9*t18*8.0_r64 &
              -nr*nzp*t5*t9*t12*t18*4.0_r64-nrp*nz*t6*t9*t12*t18*4.0_r64-nrp*nz*t4*t9*t15*t18*4.0_r64 &
              -nr*nzp*t5*t9*t15*t18*4.0_r64
              et38 = nr*nzp*t2*t6*t18*t27+nrp*nz*t4*t12*t18*t27+nr*nrp*rho*t2*t5*zh+nr*nrp*rho*t2*t7*zh*4.0_r64 &
              -nr*nrp*rhop*t2*t4*zh*6.5e+1_r64+nr*nrp*rhop*t2*t6*zh*2.0_r64+nr*nrp*rhop*t2*t27*zh+nz*nzp*rho*t2*t4*zh*4.0_r64 &
              +nz*nzp*rho*t2*t6*zh*2.0_r64-nz*nzp*rhop*t2*t5*zh*8.0_r64-nr*nrp*t6*t9*t13*zh*4.0_r64+nr*nrp*t4*t13*t27*zh &
              +nr*nrp*t5*t16*t27*zh+chi*nz*nzp*rho*t19*t27-chi*nrp*nz*t2*t9*t18*8.0_r64+chi*nr*nzp*t12*t18*t27 &
              +chi*nr*nzp*t15*t18*t27-chi*nr*nrp*rho*t2*t9*zh*8.0_r64+nr*nrp*rho*t2*t5*t27*zh+nr*nrp*rho*t2*t7*t27*zh &
              -nr*nrp*rhop*t2*t6*t9*zh*8.0_r64+nr*nrp*rhop*t2*t4*t27*zh
              coE(7) = (rhop*t21*t30*t32*t36*(et36+et37+et38))/8.0_r64
              coK(8) = rn*t21*t35*t38*(t37+t42+t48+t49+t51+t52+t53+t60+t62+t66+t67+t69+t70+t74+t75+t76+t77+t78 &
              -nr*nzp*rhop*t2*(3.0_r64*ic)+nr*nzp*rhop*t2*t4*(3.0_r64*ic))*(-1.0_r64/4.0_r64)
              coE(8) = (rn*t21*t28*t31*t35*(t40+t41+t54+t58+t63+t64+t68+nr*nzp*rho*t2*ic+nrp*nz*rho*t18*(5.0_r64*ic) &
              -nr*nzp*rhop*t18*(4.0_r64*ic)+nr*nrp*t12*zh*(5.0_r64*ic)+nr*nrp*t15*zh*(4.0_r64*ic)-chi*nr*nzp*rhop*t2*ic &
              -nr*nzp*rho*t2*t4*ic-nrp*nz*rho*t4*t18*ic+nr*nzp*rhop*t2*t5*ic-nr*nrp*t4*t12*zh*ic))/4.0_r64
              et39 = nr*nrp*1.0_r64/t34**4*4.0_r64-chi*nz*nzp*t20*8.0_r64-nrp*nz*rho*t19*5.0_r64+nr*nzp*rhop*t19*5.0_r64 &
              -nr*nrp*t12*t18*5.0_r64-nr*nrp*t15*t18*5.0_r64-nz*nzp*t2*t18*4.0_r64-nr*nrp*t4*1.0_r64/t34**4*6.0_r64 &
              +nr*nrp*t6*1.0_r64/t34**4*2.0_r64-nr*nrp*t9*1.0_r64/t34**4*4.0_r64-chi*nr*nzp*rho*t19*8.0_r64 &
              +chi*nrp*nz*rhop*t19*8.0_r64-chi*nr*nrp*t2*t12-chi*nr*nrp*t2*t15+chi*nr*nrp*t2*t18*1.7e+1_r64 &
              +chi*nz*nzp*t9*t20*3.2e+1_r64-nrp*nz*rho*t4*t19*3.0_r64+nrp*nz*rho*t9*t19*2.0e+1_r64+nr*nzp*rhop*t4*t19*3.0_r64 &
              -nr*nzp*rhop*t9*t19*2.0e+1_r64+nr*nrp*t2*t5*t12+nr*nrp*t2*t5*t15-nr*nrp*t2*t5*t18-nr*nrp*t4*t12*t18*3.0_r64 &
              -nr*nrp*t4*t15*t18*3.0_r64
              et40 = nr*nrp*t9*t12*t18*2.0e+1_r64+nr*nrp*t9*t15*t18*2.0e+1_r64+nz*nzp*t2*t4*t18*4.0_r64 &
              +nz*nzp*t2*t9*t18*1.6e+1_r64-nr*nzp*rho*t2*zh*2.0_r64+nrp*nz*rhop*t2*zh*2.0_r64+nr*nrp*t6*t27*1.0_r64/t34**4 &
              -chi*nrp*nz*rho*t2*zh*2.0_r64+chi*nr*nzp*rhop*t2*zh*2.0_r64+nrp*nz*rho*t4*t9*t19*1.2e+1_r64 &
              -nr*nzp*rhop*t4*t9*t19*1.2e+1_r64-nr*nrp*t2*t5*t9*t12*4.0_r64-nr*nrp*t2*t5*t9*t15*4.0_r64 &
              +nr*nrp*t4*t9*t12*t18*1.2e+1_r64+nr*nrp*t4*t9*t15*t18*1.2e+1_r64+nr*nrp*t2*t5*t18*t27 &
              -nz*nzp*t2*t4*t9*t18*1.6e+1_r64+nr*nzp*rho*t2*t4*zh*2.0_r64+nrp*nz*rho*t2*t5*zh*2.0_r64 &
              +nr*nzp*rho*t2*t9*zh*8.0_r64-nrp*nz*rhop*t2*t4*zh*2.0_r64-nr*nzp*rhop*t2*t5*zh*2.0_r64 &
              -nrp*nz*rhop*t2*t9*zh*8.0_r64
              et41 = chi*nr*nzp*rho*t9*t19*3.2e+1_r64-chi*nrp*nz*rhop*t9*t19*3.2e+1_r64-chi*nr*nrp*t2*t9*t18*6.8e+1_r64 &
              +chi*nr*nrp*t2*t12*t27+chi*nr*nrp*t2*t15*t27+chi*nrp*nz*rho*t2*t9*zh*8.0_r64-chi*nr*nzp*rhop*t2*t9*zh*8.0_r64 &
              -nr*nzp*rho*t2*t4*t9*zh*8.0_r64-nrp*nz*rho*t2*t5*t9*zh*8.0_r64+nrp*nz*rhop*t2*t4*t9*zh*8.0_r64 &
              +nr*nzp*rhop*t2*t5*t9*zh*8.0_r64
              coK(9) = (t21*t35*t39*(et39+et40+et41)*(-1.0_r64/8.0_r64))/rho
              et42 = -t47+nz*nzp*t20*9.0_r64+nr*nzp*rho*t19*9.0_r64-nrp*nz*rhop*t19*9.0_r64+nr*nrp*t2*t12*3.0_r64 &
              +nr*nrp*t2*t15*3.0_r64-nr*nrp*t2*t18*2.1e+1_r64+nz*nzp*t4*t20*2.3e+1_r64-nz*nzp*t9*t20*4.0_r64 &
              -chi*nr*nrp*1.0_r64/t34**4*1.0e+1_r64+nr*nrp*t5*1.0_r64/t34**4*1.2e+1_r64-nr*nrp*t7*1.0_r64/t34**4*2.0_r64 &
              +nz*nzp*t4*1.0_r64/t34**4*8.0_r64-nz*nzp*t6*1.0_r64/t34**4*4.0_r64+chi*nrp*nz*rho*t19*2.9e+1_r64 &
              -chi*nr*nzp*rhop*t19*2.9e+1_r64+chi*nr*nrp*t12*t18*2.9e+1_r64+chi*nr*nrp*t15*t18*2.9e+1_r64 &
              +chi*nz*nzp*t2*t18*1.6e+1_r64+nr*nzp*rho*t4*t19*2.3e+1_r64+nrp*nz*rho*t5*t19*3.0_r64-nr*nzp*rho*t9*t19*4.0_r64 &
              -nrp*nz*rhop*t4*t19*2.3e+1_r64-nr*nzp*rhop*t5*t19*3.0_r64+nrp*nz*rhop*t19*t27
              et43 = nr*nrp*t2*t4*t12*-2.0_r64-nr*nrp*t2*t6*t12-nr*nrp*t2*t4*t15*2.0_r64-nr*nrp*t2*t6*t15 &
              -nr*nrp*t2*t4*t18*4.4e+1_r64+nr*nrp*t2*t6*t18+nr*nrp*t5*t12*t18*3.0_r64+nr*nrp*t5*t15*t18*3.0_r64 &
              +nr*nrp*t2*t18*t27-nz*nzp*t2*t5*t18*1.6e+1_r64+nz*nzp*t4*t20*t27+nrp*nz*rho*t2*zh*6.0_r64 &
              -nr*nzp*rhop*t2*zh*6.0_r64+chi*nr*nzp*rho*t2*zh*8.0_r64-chi*nrp*nz*rhop*t2*zh*8.0_r64+nr*nzp*rho*t4*t19*t27 &
              +nrp*nz*rho*t5*t19*t27-nrp*nz*rhop*t4*t9*t19*4.0_r64-nr*nzp*rhop*t5*t9*t19*4.0_r64-nr*nrp*t2*t6*t9*t18*4.0_r64 &
              +nr*nrp*t5*t12*t18*t27+nr*nrp*t5*t15*t18*t27-nrp*nz*rho*t2*t4*zh*4.0_r64-nr*nzp*rho*t2*t5*zh*8.0_r64 &
              -nrp*nz*rho*t2*t6*zh*2.0_r64
              et44 = nr*nzp*rhop*t2*t4*zh*4.0_r64+nrp*nz*rhop*t2*t5*zh*8.0_r64+nr*nzp*rhop*t2*t6*zh*2.0_r64 &
              -chi*nrp*nz*rho*t9*t19*4.0_r64+chi*nr*nzp*rhop*t19*t27-chi*nr*nrp*t9*t12*t18*4.0_r64 &
              -chi*nr*nrp*t9*t15*t18*4.0_r64
              coE(9) = rhop*t21*t30*t32*t36*(et42+et43+et44)*(-1.0_r64/8.0_r64)
            do e = 1, 9
              ar = (e-1)/3 + 1; b = mod(e-1,3) + 1
              A((ar-1)*nt+i, (b-1)*nso+jj, md+1) = (coK(e)*vk + coE(e)*ve)*ws
            end do
            end block
          end do
        end do
      end do
      if (nk == 0) cycle
      do c = 1, npa
        if (joa(c) /= qo) cycle
        denom = tpan(qo+1) - tpan(qo)
        do i = 1, p
          tgi = tin(c) + (1.0_r64+tglp(i))/2.0_r64*(tin(c+1)-tin(c))
          rk(i) = (2.0_r64*tgi - (tpan(qo)+tpan(qo+1)))/denom
        end do
        call lagrange_interp_r64(p, tglp, p, rk, Lc)
        Ypb = matmul(Lc, xc(:,qo))
        rlo = (2.0_r64*tin(c)   - (tpan(qo)+tpan(qo+1)))/denom
        rhi = (2.0_r64*tin(c+1) - (tpan(qo)+tpan(qo+1)))/denom
        re2(1) = rlo; re2(2) = rhi
        call lagrange_interp_r64(p, tglp, 2_8, re2, IPe2)
        ec2 = matmul(IPe2, xc(:,qo)); zac = ec2(1); zbc = ec2(2)
        IPqc = matmul(IP2, Lc)
        Yq = matmul(IP2, Ypb); dYq = matmul(Dq, Yq)
        do iq = 1, qq
          spd = abs(dYq(iq)); nvq(iq) = -ic*dYq(iq)/spd; wsq(iq) = wglq(iq)*spd; wxpq(iq) = dYq(iq)*wglq(iq)
        end do
        dYp = matmul(Dp, Ypb)
        do jp = 1, p
          wsp(jp) = wglp(jp)*abs(dYp(jp))
        end do
        sumwsc = sum(wsp)
        nkc = 0
        do i = 1, nk
          ikc(i) = (abs(znear(i)-zac) + abs(znear(i)-zbc)) < 1.5_r64*sumwsc
          if (ikc(i)) then; nkc = nkc + 1; ci(nkc) = i; zc(nkc) = znear(i); zcn(nkc) = znearn(i); end if
        end do
        if (nkc > 0) then
          call axissymstok_dlpn_coef_nmode_r64(nkc, zc(1:nkc), zcn(1:nkc), qq, Yq, nvq, M, mu, &
               C1(1:3*nkc,1:3*qq,:), C2(1:3*nkc,1:3*qq,:), C3(1:3*nkc,1:3*qq,:), C4(1:3*nkc,1:3*qq,:), C5(1:3*nkc,1:3*qq,:))
          call sdspecialquad_r64(nkc, zc(1:nkc), qq, Yq, nvq, wxpq, zac, zbc, iside, &
               As(1:qq,1:nkc), Ad(1:qq,1:nkc), A1(1:qq,1:nkc), A2(1:qq,1:nkc), A3(1:qq,1:nkc), A4(1:qq,1:nkc))
          do md = 0, M
            do e = 1, 9
              ar = (e-1)/3 + 1; b = mod(e-1,3) + 1
              do iq = 1, qq
                do ia = 1, nkc
                  cc1 = C1(3*(ia-1)+ar, 3*(iq-1)+b, md+1); cc2 = C2(3*(ia-1)+ar, 3*(iq-1)+b, md+1)
                  cc3 = C3(3*(ia-1)+ar, 3*(iq-1)+b, md+1); cc4 = C4(3*(ia-1)+ar, 3*(iq-1)+b, md+1)
                  cc5 = C5(3*(ia-1)+ar, 3*(iq-1)+b, md+1)
                  Slog = cmplx(As(iq,ia),0.0_r64,r64); Dval = Ad(iq,ia)
                  Dz = cmplx(A1(iq,ia),-A2(iq,ia),r64); Dzz = cmplx(A3(iq,ia),-A4(iq,ia),r64)
                  if (mod(e,2) == 0) then
                    Acau =  ic*twopi*aimag(Dval*(cc3/nvq(iq)))
                    Ahy  = -ic*twopi*aimag(Dz  *(cc4/nvq(iq)))
                    Asup =  ic*twopi*aimag(Dzz *(cc5/nvq(iq)))
                  else
                    Acau = cmplx( twopi*real(Dval*(cc3/nvq(iq)), r64), 0.0_r64, r64)
                    Ahy  = cmplx(-twopi*real(Dz  *(cc4/nvq(iq)), r64), 0.0_r64, r64)
                    Asup = cmplx( twopi*real(Dzz *(cc5/nvq(iq)), r64), 0.0_r64, r64)
                  end if
                  Gcq(ia,iq) = twopi*cc1*Slog + cc2*wsq(iq) + Acau + Ahy + Asup
                end do
              end do
              Gpc(1:nkc,1:p) = matmul(Gcq(1:nkc,1:qq), IPqc)
              do ia = 1, nkc
                do l = 1, p
                  A((ar-1)*nt+nidx(ci(ia)), (b-1)*nso+cols+l, md+1) = &
                    A((ar-1)*nt+nidx(ci(ia)), (b-1)*nso+cols+l, md+1) + Gpc(ia,l)
                end do
              end do
            end do
          end do
        end if
        do i = 1, nk
          if (ikc(i)) cycle
          rt = real(znear(i),r64); zti = aimag(znear(i)); nrt = real(znearn(i),r64); nzt = aimag(znearn(i))
          do md = 0, M
            rn = real(md,r64)
            do jp = 1, p
              rhop = real(Ypb(jp),r64); zh = zti - aimag(Ypb(jp)); rho = rt
              nrp = real(-ic*dYp(jp)/abs(dYp(jp)),r64); nzp = aimag(-ic*dYp(jp)/abs(dYp(jp)))
              rr2 = (rho-rhop)**2 + zh*zh; chi = 1.0_r64 + rr2/(2.0_r64*rho*rhop)
              call carrier_r64(chi, md, vk, ve, Fn, An, dFn)
              block
                complex(r64) :: t2, t3, t4, t5, t7, t9, t10, t12
                complex(r64) :: t13, t15, t16, t18, t19, t21, t22, t6
                complex(r64) :: t8, t11, t14, t17, t20, t23, t27, t28
                complex(r64) :: t29, t31, t32, t34, t37, t40, t41, t42
                complex(r64) :: t43, t44, t51, t52, t53, t54, t55, t56
                complex(r64) :: t57, t58, t60, t61, t65, t66, t67, t68
                complex(r64) :: t71, t72, t73, t74, t30, t33, t35, t36
                complex(r64) :: t38, t45, t46, t47, t48, t49, t62, t63
                complex(r64) :: t64, t69, t70, t75, t76, t77, t78, t39
                complex(r64) :: et1, et2, et3, et4, t50, et5, et6, et7
                complex(r64) :: et8, et9, et10, et11, et12, et13, et14, et15
                complex(r64) :: et16, et17, et18, et19, et20, et21, et22, et23
                complex(r64) :: et24, et25, et26, et27, et28, et29, et30, et31
                complex(r64) :: et32, et33, et34, et35, et36, et37, et38, et39
                complex(r64) :: et40, et41, et42, et43, et44
                complex(r64) :: coK(9), coE(9)
                real(r64) :: nr, nz
                nr = nrt; nz = nzt
                t2 = rho*rhop
                t3 = chi+1.0_r64
                t4 = chi**2
                t5 = chi**3
                t7 = chi**5
                t9 = rn**2
                t10 = rn**3
                t12 = rho**2
                t13 = rho**3
                t15 = rhop**2
                t16 = rhop**3
                t18 = zh**2
                t19 = zh**3
                t21 = 1.0_r64/pi
                t22 = chi-1.0_r64
                t6 = t4**2
                t8 = t4**3
                t11 = t9**2
                t14 = t12**2
                t17 = t15**2
                t20 = t18**2
                t23 = 1.0_r64/t12
                t27 = t9*4.0_r64
                t28 = 1.0_r64/t3
                t29 = t4-1.0_r64
                t31 = 1.0_r64/t22**2
                t32 = 1.0_r64/t22**3
                t34 = 1.0_r64/sqrt(t2)
                t37 = nz*nzp*t19*ic
                t40 = nz*nzp*t2*zh*(2.0_r64*ic)
                t41 = chi*nz*nzp*t19*(4.0_r64*ic)
                t42 = nr*nzp*rho*t18*ic
                t43 = nrp*nz*rhop*t18*ic
                t44 = nr*nrp*t2*zh*(1.0e+1_r64*ic)
                t51 = chi*nr*nrp*t12*zh*ic
                t52 = chi*nr*nrp*t15*zh*ic
                t53 = chi*nrp*nz*rho*t18*ic
                t54 = chi*nr*nzp*rho*t18*(4.0_r64*ic)
                t55 = chi*nr*nzp*rhop*t18*ic
                t56 = chi*nrp*nz*rhop*t18*(4.0_r64*ic)
                t57 = chi*nr*nrp*t2*zh*(1.2e+1_r64*ic)
                t58 = nr*nrp*t2*t5*zh*(4.0_r64*ic)
                t60 = nr*nrp*t2*t9*zh*(4.0_r64*ic)
                t61 = nz*nzp*t9*t19*(4.0_r64*ic)
                t65 = nr*nzp*rho*t9*t18*(4.0_r64*ic)
                t66 = nrp*nz*rhop*t9*t18*(4.0_r64*ic)
                t67 = nr*nrp*t2*t4*zh*(8.0_r64*ic)
                t68 = nz*nzp*t2*t4*zh*-(2.0_r64*ic)
                t71 = chi*nr*nrp*t9*t12*zh*(4.0_r64*ic)
                t72 = chi*nr*nrp*t9*t15*zh*(4.0_r64*ic)
                t73 = chi*nrp*nz*rho*t9*t18*(4.0_r64*ic)
                t74 = chi*nr*nzp*rhop*t9*t18*(4.0_r64*ic)
                t30 = t28**2
                t33 = t27-1.0_r64
                t35 = t34**5
                t36 = t34**7
                t38 = 1.0_r64/t29
                t45 = nr*nzp*1.0_r64/t34**4*4.0_r64
                t46 = nrp*nz*1.0_r64/t34**4*4.0_r64
                t47 = nz*nzp*1.0_r64/t34**4*4.0_r64
                t48 = -t43
                t49 = -t44
                t62 = -t55
                t63 = -t56
                t64 = -t57
                t69 = -t61
                t70 = -t65
                t75 = t4*t60
                t76 = -t71
                t77 = -t72
                t78 = -t73
                t39 = t38**2
                et1 = t47-chi*nr*nrp*t14*1.0e+1_r64-chi*nr*nrp*t17*1.0e+1_r64+nr*nrp*t2*t12*2.8e+1_r64+nr*nrp*t2*t15*2.8e+1_r64 &
                +nr*nrp*t5*t14*2.0_r64+nr*nrp*t5*t17*2.0_r64-nz*nzp*t12*t18*5.0_r64-nz*nzp*t15*t18*5.0_r64-nr*nzp*t13*zh*5.0_r64 &
                +nrp*nz*t16*zh*5.0_r64-chi*nr*nrp*1.0_r64/t34**4*6.8e+1_r64+nr*nrp*t5*1.0_r64/t34**4*3.6e+1_r64 &
                -nr*nrp*t7*1.0_r64/t34**4*1.6e+1_r64-nz*nzp*t4*1.0_r64/t34**4*6.0_r64+nz*nzp*t6*1.0_r64/t34**4*2.0_r64 &
                -nz*nzp*t9*1.0_r64/t34**4*4.0_r64+chi*nr*nrp*t9*t14*4.0e+1_r64+chi*nr*nrp*t9*t17*4.0e+1_r64-chi*nz*nzp*t2*t12 &
                -chi*nz*nzp*t2*t15+chi*nz*nzp*t2*t18*1.7e+1_r64-chi*nrp*nz*t13*zh*1.0e+1_r64+chi*nr*nzp*t16*zh*1.0e+1_r64 &
                -nr*nrp*t2*t4*t12*4.0_r64
                et2 = nr*nrp*t2*t6*t12*8.0_r64-nr*nrp*t2*t4*t15*4.0_r64+nr*nrp*t2*t6*t15*8.0_r64-nr*nrp*t2*t9*t12*4.0e+1_r64 &
                -nr*nrp*t2*t9*t15*4.0e+1_r64-nr*nrp*t5*t9*t14*8.0_r64-nr*nrp*t5*t9*t17*8.0_r64+nz*nzp*t2*t5*t12+nz*nzp*t2*t5*t15 &
                -nz*nzp*t2*t5*t18-nz*nzp*t4*t12*t18*3.0_r64-nz*nzp*t4*t15*t18*3.0_r64+nz*nzp*t9*t12*t18*2.0e+1_r64 &
                +nz*nzp*t9*t15*t18*2.0e+1_r64+nrp*nz*rho*t2*zh*1.9e+1_r64-nr*nzp*rhop*t2*zh*1.9e+1_r64-nr*nzp*t4*t13*zh*3.0_r64 &
                +nrp*nz*t5*t13*zh*2.0_r64+nrp*nz*t4*t16*zh*3.0_r64-nr*nzp*t5*t16*zh*2.0_r64+nr*nzp*t9*t13*zh*2.0e+1_r64 &
                -nrp*nz*t9*t16*zh*2.0e+1_r64+chi*nr*nrp*t9*1.0_r64/t34**4*1.76e+2_r64+nr*nrp*t5*t9*1.0_r64/t34**4*4.8e+1_r64
                et3 = nr*nrp*t7*t9*1.0_r64/t34**4*-3.2e+1_r64+nz*nzp*t6*t27*1.0_r64/t34**4+chi*nr*nzp*rho*t2*zh*2.6e+1_r64 &
                -chi*nrp*nz*rhop*t2*zh*2.6e+1_r64+chi*nrp*nz*t9*t13*zh*4.0e+1_r64-chi*nr*nzp*t9*t16*zh*4.0e+1_r64 &
                -nr*nrp*t2*t4*t9*t12*1.28e+2_r64+nr*nrp*t2*t6*t9*t12*4.0e+1_r64-nr*nrp*t2*t4*t9*t15*1.28e+2_r64 &
                +nr*nrp*t2*t6*t9*t15*4.0e+1_r64-nz*nzp*t2*t5*t9*t12*4.0_r64-nz*nzp*t2*t5*t9*t15*4.0_r64 &
                +nz*nzp*t4*t9*t12*t18*1.2e+1_r64+nz*nzp*t4*t9*t15*t18*1.2e+1_r64+nz*nzp*t2*t5*t18*t27+nrp*nz*rho*t2*t4*zh &
                -nr*nzp*rho*t2*t5*zh*2.0_r64+nrp*nz*rho*t2*t6*zh*4.0_r64-nrp*nz*rho*t2*t9*zh*4.0e+1_r64-nr*nzp*rhop*t2*t4*zh &
                +nrp*nz*rhop*t2*t5*zh*2.0_r64
                et4 = nr*nzp*rhop*t2*t6*zh*-4.0_r64+nr*nzp*rhop*t2*t9*zh*4.0e+1_r64+nr*nzp*t4*t9*t13*zh*1.2e+1_r64 &
                -nrp*nz*t5*t9*t13*zh*8.0_r64-nrp*nz*t4*t9*t16*zh*1.2e+1_r64+nr*nzp*t5*t9*t16*zh*8.0_r64 &
                -chi*nz*nzp*t2*t9*t18*6.8e+1_r64+chi*nz*nzp*t2*t12*t27+chi*nz*nzp*t2*t15*t27-chi*nr*nzp*rho*t2*t9*zh*1.04e+2_r64 &
                +chi*nrp*nz*rhop*t2*t9*zh*1.04e+2_r64-nrp*nz*rho*t2*t4*t9*zh*7.6e+1_r64+nr*nzp*rho*t2*t5*t9*zh*8.0_r64 &
                +nrp*nz*rho*t2*t6*t9*zh*2.0e+1_r64+nr*nzp*rhop*t2*t4*t9*zh*7.6e+1_r64-nrp*nz*rhop*t2*t5*t9*zh*8.0_r64 &
                -nr*nzp*rhop*t2*t6*t9*zh*2.0e+1_r64
                coK(1) = rhop*t21*t36*t39*(et1+et2+et3+et4)*(-1.0_r64/8.0_r64)
                t50 = 1.0_r64/t33
                et5 = nr*nrp*1.0_r64/t34**4*5.0e+1_r64+nr*nrp*t14*1.5e+1_r64+nr*nrp*t17*1.5e+1_r64+nr*nrp*t4*t14*1.9e+1_r64 &
                -nr*nrp*t6*t14*2.0_r64+nr*nrp*t4*t17*1.9e+1_r64-nr*nrp*t6*t17*2.0_r64-nr*nrp*t9*t14*6.0e+1_r64 &
                -nr*nrp*t9*t17*6.0e+1_r64+nz*nzp*t2*t12*3.0_r64+nz*nzp*t2*t15*3.0_r64-nz*nzp*t2*t18*2.1e+1_r64 &
                +nrp*nz*t13*zh*1.5e+1_r64-nr*nzp*t16*zh*1.5e+1_r64-chi*nz*nzp*1.0_r64/t34**4*1.0e+1_r64 &
                +nr*nrp*t4*1.0_r64/t34**4*1.78e+2_r64-nr*nrp*t6*1.0_r64/t34**4*5.2e+1_r64+nr*nrp*t8*1.0_r64/t34**4*1.6e+1_r64 &
                -nr*nrp*t9*1.0_r64/t34**4*2.16e+2_r64+nr*nrp*t11*1.0_r64/t34**4*1.6e+1_r64+nz*nzp*t5*1.0_r64/t34**4*1.2e+1_r64 &
                -nz*nzp*t7*1.0_r64/t34**4*2.0_r64-chi*nr*nrp*t2*t12*1.24e+2_r64
                et6 = chi*nr*nrp*t2*t15*-1.24e+2_r64+chi*nz*nzp*t12*t18*2.9e+1_r64+chi*nz*nzp*t15*t18*2.9e+1_r64 &
                +chi*nr*nzp*t13*zh*2.9e+1_r64-chi*nrp*nz*t16*zh*2.9e+1_r64+nr*nrp*t2*t5*t12*4.0_r64-nr*nrp*t2*t7*t12*8.0_r64 &
                +nr*nrp*t2*t5*t15*4.0_r64-nr*nrp*t2*t7*t15*8.0_r64-nr*nrp*t4*t9*t14*8.0e+1_r64+nr*nrp*t4*t11*t14*1.6e+1_r64 &
                +nr*nrp*t6*t9*t14*1.2e+1_r64-nr*nrp*t4*t9*t17*8.0e+1_r64-nr*nrp*t6*t11*t14*1.6e+1_r64 &
                +nr*nrp*t4*t11*t17*1.6e+1_r64+nr*nrp*t6*t9*t17*1.2e+1_r64-nr*nrp*t6*t11*t17*1.6e+1_r64-nz*nzp*t2*t4*t12*2.0_r64 &
                -nz*nzp*t2*t6*t12-nz*nzp*t2*t4*t15*2.0_r64-nz*nzp*t2*t6*t15-nz*nzp*t2*t9*t12*1.2e+1_r64 &
                -nz*nzp*t2*t4*t18*4.4e+1_r64
                et7 = nz*nzp*t2*t6*t18-nz*nzp*t2*t9*t15*1.2e+1_r64+nz*nzp*t2*t9*t18*8.8e+1_r64-nz*nzp*t2*t11*t18*1.6e+1_r64 &
                +nz*nzp*t5*t12*t18*3.0_r64+nz*nzp*t5*t15*t18*3.0_r64-nr*nzp*rho*t2*zh*3.3e+1_r64+nrp*nz*rhop*t2*zh*3.3e+1_r64 &
                +nrp*nz*t4*t13*zh*1.9e+1_r64+nr*nzp*t5*t13*zh*3.0_r64-nrp*nz*t6*t13*zh*2.0_r64-nr*nzp*t4*t16*zh*1.9e+1_r64 &
                -nrp*nz*t5*t16*zh*3.0_r64+nr*nzp*t6*t16*zh*2.0_r64-nrp*nz*t9*t13*zh*6.0e+1_r64+nr*nzp*t9*t16*zh*6.0e+1_r64 &
                +chi*nz*nzp*t9*1.0_r64/t34**4*4.0e+1_r64-nr*nrp*t4*t9*1.0_r64/t34**4*6.88e+2_r64 &
                +nr*nrp*t4*t11*1.0_r64/t34**4*4.8e+1_r64+nr*nrp*t6*t9*1.0_r64/t34**4*1.84e+2_r64 &
                -nr*nrp*t6*t11*1.0_r64/t34**4*4.8e+1_r64-nr*nrp*t8*t9*1.0_r64/t34**4*4.8e+1_r64
                et8 = nr*nrp*t8*t11*1.0_r64/t34**4*-1.6e+1_r64-nz*nzp*t5*t9*1.0_r64/t34**4*4.8e+1_r64 &
                +nz*nzp*t7*t9*1.0_r64/t34**4*8.0_r64-chi*nrp*nz*rho*t2*zh*9.1e+1_r64+chi*nr*nzp*rhop*t2*zh*9.1e+1_r64 &
                -chi*nr*nzp*t9*t13*zh*1.2e+2_r64+chi*nr*nzp*t11*t13*zh*1.6e+1_r64+chi*nrp*nz*t9*t16*zh*1.2e+2_r64 &
                -chi*nrp*nz*t11*t16*zh*1.6e+1_r64-nr*nrp*t2*t5*t9*t12*1.6e+1_r64+nr*nrp*t2*t7*t9*t12*2.4e+1_r64 &
                -nr*nrp*t2*t5*t9*t15*1.6e+1_r64+nr*nrp*t2*t7*t11*t12*3.2e+1_r64+nr*nrp*t2*t7*t9*t15*2.4e+1_r64 &
                +nr*nrp*t2*t7*t11*t15*3.2e+1_r64+nz*nzp*t2*t4*t9*t12*8.0_r64+nz*nzp*t2*t4*t9*t15*8.0_r64 &
                +nz*nzp*t2*t4*t9*t18*1.76e+2_r64-nz*nzp*t2*t6*t9*t18*8.0_r64+nz*nzp*t2*t6*t11*t18*1.6e+1_r64
                et9 = nz*nzp*t5*t9*t12*t18*-8.0_r64-nz*nzp*t5*t11*t12*t18*1.6e+1_r64+nz*nzp*t2*t6*t12*t27 &
                -nz*nzp*t5*t9*t15*t18*8.0_r64-nz*nzp*t5*t11*t15*t18*1.6e+1_r64+nz*nzp*t2*t6*t15*t27 &
                -nr*nzp*rho*t2*t4*zh*6.5e+1_r64-nrp*nz*rho*t2*t5*zh+nr*nzp*rho*t2*t6*zh*2.0_r64-nrp*nz*rho*t2*t7*zh*4.0_r64 &
                +nr*nzp*rho*t2*t9*zh*1.36e+2_r64-nr*nzp*rho*t2*t11*zh*1.6e+1_r64+nrp*nz*rhop*t2*t4*zh*6.5e+1_r64 &
                +nr*nzp*rhop*t2*t5*zh-nrp*nz*rhop*t2*t6*zh*2.0_r64+nr*nzp*rhop*t2*t7*zh*4.0_r64-nrp*nz*rhop*t2*t9*zh*1.36e+2_r64 &
                +nrp*nz*rhop*t2*t11*zh*1.6e+1_r64-nrp*nz*t4*t9*t13*zh*8.0e+1_r64-nr*nzp*t5*t9*t13*zh*8.0_r64 &
                +nrp*nz*t4*t11*t13*zh*1.6e+1_r64
                et10 = nrp*nz*t6*t9*t13*zh*1.2e+1_r64+nr*nzp*t4*t9*t16*zh*8.0e+1_r64-nr*nzp*t5*t11*t13*zh*1.6e+1_r64 &
                +nrp*nz*t5*t9*t16*zh*8.0_r64-nrp*nz*t6*t11*t13*zh*1.6e+1_r64-nr*nzp*t4*t11*t16*zh*1.6e+1_r64 &
                -nr*nzp*t6*t9*t16*zh*1.2e+1_r64+nrp*nz*t5*t11*t16*zh*1.6e+1_r64+nr*nzp*t6*t11*t16*zh*1.6e+1_r64 &
                +chi*nr*nrp*t2*t9*t12*5.04e+2_r64-chi*nr*nrp*t2*t11*t12*3.2e+1_r64+chi*nr*nrp*t2*t9*t15*5.04e+2_r64 &
                -chi*nr*nrp*t2*t11*t15*3.2e+1_r64-chi*nz*nzp*t9*t12*t18*1.2e+2_r64+chi*nz*nzp*t11*t12*t18*1.6e+1_r64 &
                -chi*nz*nzp*t9*t15*t18*1.2e+2_r64+chi*nz*nzp*t11*t15*t18*1.6e+1_r64+chi*nrp*nz*rho*t2*t9*zh*3.72e+2_r64 &
                -chi*nrp*nz*rho*t2*t11*zh*3.2e+1_r64-chi*nr*nzp*rhop*t2*t9*zh*3.72e+2_r64
                et11 = chi*nr*nzp*rhop*t2*t11*zh*3.2e+1_r64+nr*nzp*rho*t2*t4*t9*zh*2.64e+2_r64-nr*nzp*rho*t2*t4*t11*zh*1.6e+1_r64 &
                -nr*nzp*rho*t2*t6*t9*zh*1.6e+1_r64+nrp*nz*rho*t2*t5*t11*zh*1.6e+1_r64+nrp*nz*rho*t2*t7*t9*zh*1.2e+1_r64 &
                +nr*nzp*rho*t2*t6*t11*zh*3.2e+1_r64+nrp*nz*rho*t2*t7*t11*zh*1.6e+1_r64-nrp*nz*rhop*t2*t4*t9*zh*2.64e+2_r64 &
                +nrp*nz*rhop*t2*t4*t11*zh*1.6e+1_r64+nrp*nz*rhop*t2*t6*t9*zh*1.6e+1_r64-nr*nzp*rhop*t2*t5*t11*zh*1.6e+1_r64 &
                -nr*nzp*rhop*t2*t7*t9*zh*1.2e+1_r64-nrp*nz*rhop*t2*t6*t11*zh*3.2e+1_r64-nr*nzp*rhop*t2*t7*t11*zh*1.6e+1_r64
                coE(1) = (t21*t30*t32*t35*t50*(et5+et6+et7+et8+et9+et10+et11))/(rho*8.0_r64)
                et12 = nr*nrp*t16*-(9.0_r64*ic)+chi*nr*nrp*t13*ic-nr*nrp*rho*t2*(1.9e+1_r64*ic)+nz*nzp*rho*t18*ic &
                -nz*nzp*rhop*t2*(3.0_r64*ic)+nr*nrp*t4*t16*(8.0_r64*ic)-nrp*nz*t2*zh*(1.0e+1_r64*ic)+nr*nzp*t12*zh*ic &
                +nr*nzp*t15*zh*(9.0_r64*ic)+chi*nr*nrp*rhop*t2*(2.3e+1_r64*ic)-chi*nz*nzp*rhop*t18*ic &
                -chi*nr*nrp*t9*t13*(4.0_r64*ic)-chi*nr*nzp*t2*zh*(2.0_r64*ic)+chi*nrp*nz*t12*zh*ic+chi*nrp*nz*t15*zh*ic &
                +nr*nrp*rho*t2*t4*(1.6e+1_r64*ic)+nr*nrp*rho*t2*t9*(4.0_r64*ic)-nr*nrp*rhop*t2*t5*(2.0e+1_r64*ic) &
                -nz*nzp*rho*t9*t18*(4.0_r64*ic)+nz*nzp*rhop*t2*t4*(3.0_r64*ic)+nr*nrp*t4*t9*t16*(4.0_r64*ic)
                et13 = nrp*nz*t2*t4*zh*(8.0_r64*ic)+nrp*nz*t2*t9*zh*(4.0_r64*ic)-nr*nzp*t4*t15*zh*(8.0_r64*ic) &
                -nr*nzp*t9*t12*zh*(4.0_r64*ic)+chi*nr*nzp*t2*t9*zh*(8.0_r64*ic)-chi*nrp*nz*t9*t12*zh*(4.0_r64*ic) &
                -chi*nrp*nz*t9*t15*zh*(4.0_r64*ic)+nr*nrp*rho*t2*t4*t9*(8.0_r64*ic)-nr*nrp*rhop*t2*t5*t9*(4.0_r64*ic) &
                +nrp*nz*t2*t4*t9*zh*(4.0_r64*ic)-nr*nzp*t4*t9*t15*zh*(4.0_r64*ic)-chi*nr*nrp*rhop*t2*t9*(8.0_r64*ic) &
                +chi*nz*nzp*rhop*t9*t18*(4.0_r64*ic)
                coK(2) = rn*t21*t35*t38*(et12+et13)*(-1.0_r64/4.0_r64)
                et14 = rn*nr*nrp*t13*(5.0_r64*ic)-nr*nrp*t10*t13*(2.0e+1_r64*ic)-chi*rn*nr*nrp*t16*(8.0_r64*ic) &
                +chi*nr*nrp*t10*t16*(3.2e+1_r64*ic)+rn*nr*nrp*rhop*t2*(1.1e+1_r64*ic)+rn*nz*nzp*rho*t2*ic &
                -rn*nz*nzp*rhop*t18*(4.0_r64*ic)-rn*nr*nrp*t4*t13*ic+rn*nr*nrp*t5*t16*(4.0_r64*ic) &
                -nr*nrp*rhop*t2*t10*(3.2e+1_r64*ic)-nz*nzp*rho*t2*t10*(4.0_r64*ic)+nz*nzp*rhop*t10*t18*(1.6e+1_r64*ic) &
                +nr*nrp*t4*t10*t13*(4.0_r64*ic)-nr*nrp*t5*t10*t16*(1.6e+1_r64*ic)-rn*nr*nzp*t2*zh*(8.0_r64*ic) &
                +rn*nrp*nz*t12*zh*(5.0_r64*ic)+rn*nrp*nz*t15*zh*(4.0_r64*ic)+nr*nzp*t2*t10*zh*(3.2e+1_r64*ic) &
                -nrp*nz*t10*t12*zh*(2.0e+1_r64*ic)-nrp*nz*t10*t15*zh*(1.6e+1_r64*ic)
                et15 = chi*nrp*nz*t2*t10*zh*(4.8e+1_r64*ic)-chi*nr*nzp*t10*t12*zh*(1.6e+1_r64*ic) &
                -chi*nr*nzp*t10*t15*zh*(3.2e+1_r64*ic)+rn*nr*nrp*rho*t2*t5*(8.0_r64*ic)+rn*nr*nrp*rhop*t2*t4*(5.0_r64*ic) &
                -rn*nr*nrp*rhop*t2*t6*(4.0_r64*ic)-rn*nz*nzp*rho*t2*t4*ic+rn*nz*nzp*rhop*t2*t5*ic &
                -nr*nrp*rho*t2*t5*t10*(3.2e+1_r64*ic)-nr*nrp*rhop*t2*t4*t10*(4.4e+1_r64*ic)+nr*nrp*rhop*t2*t6*t10*(2.8e+1_r64*ic) &
                +nz*nzp*rho*t2*t4*t10*(4.0_r64*ic)-nz*nzp*rhop*t2*t5*t10*(4.0_r64*ic)+rn*nrp*nz*t2*t5*zh*(4.0_r64*ic) &
                -rn*nrp*nz*t4*t12*zh*ic-rn*nr*nzp*t5*t15*zh*(4.0_r64*ic)-nrp*nz*t2*t5*t10*zh*(1.6e+1_r64*ic) &
                +nrp*nz*t4*t10*t12*zh*(4.0_r64*ic)
                et16 = nr*nzp*t5*t10*t15*zh*(1.6e+1_r64*ic)-chi*rn*nr*nrp*rho*t2*(2.0e+1_r64*ic) &
                +chi*rn*nz*nzp*rho*t18*(4.0_r64*ic)-chi*rn*nz*nzp*rhop*t2*ic+chi*nr*nrp*rho*t2*t10*(8.0e+1_r64*ic) &
                -chi*nz*nzp*rho*t10*t18*(1.6e+1_r64*ic)+chi*nz*nzp*rhop*t2*t10*(4.0_r64*ic)-chi*rn*nrp*nz*t2*zh*(1.2e+1_r64*ic) &
                +chi*rn*nr*nzp*t12*zh*(4.0_r64*ic)+chi*rn*nr*nzp*t15*zh*(8.0_r64*ic)
                coE(2) = t21*t28*t31*t35*t50*(et14+et15+et16)*(-1.0_r64/4.0_r64)
                et17 = -t46-nz*nzp*rhop*t19*5.0_r64-nr*nzp*t2*t18*1.0e+1_r64+nrp*nz*t12*t18*5.0_r64+nrp*nz*t15*t18*5.0_r64 &
                +nr*nrp*t13*zh*5.0_r64+nrp*nz*t4*1.0_r64/t34**4*6.0_r64-nrp*nz*t6*1.0_r64/t34**4*2.0_r64 &
                +nrp*nz*t27*1.0_r64/t34**4+chi*nz*nzp*rho*t19*8.0_r64+chi*nrp*nz*t2*t12+chi*nrp*nz*t2*t15 &
                -chi*nrp*nz*t2*t18*1.7e+1_r64+chi*nr*nzp*t12*t18*8.0_r64+chi*nr*nzp*t15*t18*1.0e+1_r64 &
                -chi*nr*nrp*t16*zh*1.0e+1_r64-nz*nzp*rhop*t4*t19*3.0_r64+nz*nzp*rhop*t9*t19*2.0e+1_r64-nrp*nz*t2*t5*t12 &
                -nrp*nz*t2*t5*t15-nr*nzp*t2*t4*t18*6.0_r64+nrp*nz*t2*t5*t18+nr*nzp*t2*t9*t18*4.0e+1_r64+nrp*nz*t4*t12*t18*3.0_r64 &
                +nrp*nz*t4*t15*t18*3.0_r64-nr*nzp*t5*t15*t18*2.0_r64
                et18 = nrp*nz*t9*t12*t18*-2.0e+1_r64-nrp*nz*t9*t15*t18*2.0e+1_r64+nr*nrp*rhop*t2*zh*1.9e+1_r64 &
                +nz*nzp*rho*t2*zh*2.0_r64+nr*nrp*t4*t13*zh*3.0_r64+nr*nrp*t5*t16*zh*2.0_r64-nr*nrp*t9*t13*zh*2.0e+1_r64 &
                -nrp*nz*t6*t9*1.0_r64/t34**4*4.0_r64-chi*nr*nrp*rho*t2*zh*2.6e+1_r64-chi*nz*nzp*rhop*t2*zh*2.0_r64 &
                +chi*nr*nrp*t9*t16*zh*4.0e+1_r64+nz*nzp*rhop*t4*t9*t19*1.2e+1_r64+nr*nzp*t2*t4*t9*t18*2.4e+1_r64 &
                -nrp*nz*t2*t5*t9*t18*4.0_r64-nrp*nz*t4*t9*t12*t18*1.2e+1_r64+nrp*nz*t2*t5*t12*t27-nrp*nz*t4*t9*t15*t18*1.2e+1_r64 &
                +nr*nzp*t5*t9*t15*t18*8.0_r64+nrp*nz*t2*t5*t15*t27+nr*nrp*rho*t2*t5*zh*2.0_r64+nr*nrp*rhop*t2*t4*zh &
                +nr*nrp*rhop*t2*t6*zh*4.0_r64-nr*nrp*rhop*t2*t9*zh*4.0e+1_r64
                et19 = nz*nzp*rho*t2*t4*zh*-2.0_r64-nz*nzp*rho*t2*t9*zh*8.0_r64+nz*nzp*rhop*t2*t5*zh*2.0_r64 &
                -nr*nrp*t4*t9*t13*zh*1.2e+1_r64-nr*nrp*t5*t9*t16*zh*8.0_r64-chi*nz*nzp*rho*t9*t19*3.2e+1_r64 &
                -chi*nrp*nz*t2*t9*t12*4.0_r64-chi*nrp*nz*t2*t9*t15*4.0_r64+chi*nrp*nz*t2*t9*t18*6.8e+1_r64 &
                -chi*nr*nzp*t9*t12*t18*3.2e+1_r64-chi*nr*nzp*t9*t15*t18*4.0e+1_r64+chi*nr*nrp*rho*t2*t9*zh*1.04e+2_r64 &
                +chi*nz*nzp*rhop*t2*t9*zh*8.0_r64-nr*nrp*rho*t2*t5*t9*zh*8.0_r64-nr*nrp*rhop*t2*t4*t9*zh*7.6e+1_r64 &
                +nr*nrp*rhop*t2*t6*t9*zh*2.0e+1_r64+nz*nzp*rho*t2*t4*t9*zh*8.0_r64-nz*nzp*rhop*t2*t5*t9*zh*8.0_r64
                coK(3) = (rhop*t21*t36*t39*(et17+et18+et19))/8.0_r64
                et20 = -t45+nz*nzp*rho*t19*9.0_r64+nrp*nz*t2*t12*3.0_r64+nrp*nz*t2*t15*3.0_r64-nrp*nz*t2*t18*2.1e+1_r64 &
                +nr*nzp*t12*t18*9.0_r64+nr*nzp*t15*t18*1.5e+1_r64-nr*nrp*t16*zh*1.5e+1_r64-chi*nrp*nz*1.0_r64/t34**4*1.0e+1_r64 &
                +nr*nzp*t4*1.0_r64/t34**4*8.0_r64+nrp*nz*t5*1.0_r64/t34**4*1.2e+1_r64-nr*nzp*t6*1.0_r64/t34**4*4.0_r64 &
                -nrp*nz*t7*1.0_r64/t34**4*2.0_r64-chi*nz*nzp*rhop*t19*2.9e+1_r64-chi*nr*nzp*t2*t18*5.8e+1_r64 &
                +chi*nrp*nz*t12*t18*2.9e+1_r64+chi*nrp*nz*t15*t18*2.9e+1_r64+chi*nr*nrp*t13*zh*2.9e+1_r64 &
                +nz*nzp*rho*t4*t19*2.3e+1_r64-nz*nzp*rho*t9*t19*4.0_r64-nz*nzp*rhop*t5*t19*3.0_r64-nrp*nz*t2*t4*t12*2.0_r64 &
                -nrp*nz*t2*t6*t12-nrp*nz*t2*t4*t15*2.0_r64
                et21 = -nrp*nz*t2*t6*t15-nrp*nz*t2*t4*t18*4.4e+1_r64-nr*nzp*t2*t5*t18*6.0_r64+nrp*nz*t2*t6*t18 &
                +nr*nzp*t4*t12*t18*2.3e+1_r64+nrp*nz*t5*t12*t18*3.0_r64+nr*nzp*t4*t15*t18*1.9e+1_r64+nrp*nz*t5*t15*t18*3.0_r64 &
                -nr*nzp*t6*t15*t18*2.0_r64-nr*nzp*t9*t12*t18*4.0_r64+nrp*nz*t2*t18*t27-nr*nrp*rho*t2*zh*3.3e+1_r64 &
                -nz*nzp*rhop*t2*zh*6.0_r64+nr*nrp*t5*t13*zh*3.0_r64-nr*nrp*t4*t16*zh*1.9e+1_r64+nr*nrp*t6*t16*zh*2.0_r64 &
                +chi*nr*nrp*rhop*t2*zh*9.1e+1_r64+chi*nz*nzp*rho*t2*zh*8.0_r64-chi*nr*nrp*t9*t13*zh*4.0_r64+nz*nzp*rho*t4*t19*t27 &
                -nz*nzp*rhop*t5*t9*t19*4.0_r64-nr*nzp*t2*t5*t9*t18*8.0_r64-nrp*nz*t2*t6*t9*t18*4.0_r64 &
                -nr*nzp*t4*t9*t15*t18*4.0_r64
                et22 = nr*nzp*t4*t12*t18*t27+nrp*nz*t5*t12*t18*t27+nrp*nz*t5*t15*t18*t27+nr*nzp*t6*t15*t18*t27 &
                -nr*nrp*rho*t2*t4*zh*6.5e+1_r64+nr*nrp*rho*t2*t6*zh*2.0_r64+nr*nrp*rho*t2*t27*zh+nr*nrp*rhop*t2*t5*zh &
                +nr*nrp*rhop*t2*t7*zh*4.0_r64-nz*nzp*rho*t2*t5*zh*8.0_r64+nz*nzp*rhop*t2*t4*zh*4.0_r64 &
                +nz*nzp*rhop*t2*t6*zh*2.0_r64-nr*nrp*t6*t9*t16*zh*4.0_r64+nr*nrp*t5*t13*t27*zh+nr*nrp*t4*t16*t27*zh &
                +chi*nz*nzp*rhop*t19*t27+chi*nr*nzp*t2*t9*t18*8.0_r64-chi*nrp*nz*t9*t12*t18*4.0_r64-chi*nrp*nz*t9*t15*t18*4.0_r64 &
                -chi*nr*nrp*rhop*t2*t9*zh*8.0_r64-nr*nrp*rho*t2*t6*t9*zh*8.0_r64+nr*nrp*rho*t2*t4*t27*zh+nr*nrp*rhop*t2*t5*t27*zh &
                +nr*nrp*rhop*t2*t7*t27*zh
                coE(3) = rhop*t21*t30*t32*t36*(et20+et21+et22)*(-1.0_r64/8.0_r64)
                et23 = nr*nrp*t13*(9.0_r64*ic)-chi*nr*nrp*t16*ic+nr*nrp*rhop*t2*(1.9e+1_r64*ic)+nz*nzp*rho*t2*(3.0_r64*ic) &
                -nz*nzp*rhop*t18*ic-nr*nrp*t4*t13*(8.0_r64*ic)-nr*nzp*t2*zh*(1.0e+1_r64*ic)+nrp*nz*t12*zh*(9.0_r64*ic) &
                +nrp*nz*t15*zh*ic-chi*nr*nrp*rho*t2*(2.3e+1_r64*ic)+chi*nz*nzp*rho*t18*ic+chi*nr*nrp*t9*t16*(4.0_r64*ic) &
                -chi*nrp*nz*t2*zh*(2.0_r64*ic)+chi*nr*nzp*t12*zh*ic+chi*nr*nzp*t15*zh*ic+nr*nrp*rho*t2*t5*(2.0e+1_r64*ic) &
                -nr*nrp*rhop*t2*t4*(1.6e+1_r64*ic)-nr*nrp*rhop*t2*t9*(4.0_r64*ic)-nz*nzp*rho*t2*t4*(3.0_r64*ic) &
                +nz*nzp*rhop*t9*t18*(4.0_r64*ic)-nr*nrp*t4*t9*t13*(4.0_r64*ic)
                et24 = nr*nzp*t2*t4*zh*(8.0_r64*ic)+nr*nzp*t2*t9*zh*(4.0_r64*ic)-nrp*nz*t4*t12*zh*(8.0_r64*ic) &
                -nrp*nz*t9*t15*zh*(4.0_r64*ic)+chi*nrp*nz*t2*t9*zh*(8.0_r64*ic)-chi*nr*nzp*t9*t12*zh*(4.0_r64*ic) &
                -chi*nr*nzp*t9*t15*zh*(4.0_r64*ic)+nr*nrp*rho*t2*t5*t9*(4.0_r64*ic)-nr*nrp*rhop*t2*t4*t9*(8.0_r64*ic) &
                +nr*nzp*t2*t4*t9*zh*(4.0_r64*ic)-nrp*nz*t4*t9*t12*zh*(4.0_r64*ic)+chi*nr*nrp*rho*t2*t9*(8.0_r64*ic) &
                -chi*nz*nzp*rho*t9*t18*(4.0_r64*ic)
                coK(4) = rn*t15*t21*t36*t38*(et23+et24)*(-1.0_r64/4.0_r64)
                et25 = rn*nr*nrp*t16*(5.0_r64*ic)-nr*nrp*t10*t16*(2.0e+1_r64*ic)-chi*rn*nr*nrp*t13*(8.0_r64*ic) &
                +chi*nr*nrp*t10*t13*(3.2e+1_r64*ic)+rn*nr*nrp*rho*t2*(1.1e+1_r64*ic)-rn*nz*nzp*rho*t18*(4.0_r64*ic) &
                +rn*nz*nzp*rhop*t2*ic+rn*nr*nrp*t5*t13*(4.0_r64*ic)-rn*nr*nrp*t4*t16*ic-nr*nrp*rho*t2*t10*(3.2e+1_r64*ic) &
                +nz*nzp*rho*t10*t18*(1.6e+1_r64*ic)-nz*nzp*rhop*t2*t10*(4.0_r64*ic)-nr*nrp*t5*t10*t13*(1.6e+1_r64*ic) &
                +nr*nrp*t4*t10*t16*(4.0_r64*ic)+rn*nrp*nz*t2*zh*(8.0_r64*ic)-rn*nr*nzp*t12*zh*(4.0_r64*ic) &
                -rn*nr*nzp*t15*zh*(5.0_r64*ic)-nrp*nz*t2*t10*zh*(3.2e+1_r64*ic)+nr*nzp*t10*t12*zh*(1.6e+1_r64*ic) &
                +nr*nzp*t10*t15*zh*(2.0e+1_r64*ic)
                et26 = chi*nr*nzp*t2*t10*zh*-(4.8e+1_r64*ic)+chi*nrp*nz*t10*t12*zh*(3.2e+1_r64*ic) &
                +chi*nrp*nz*t10*t15*zh*(1.6e+1_r64*ic)+rn*nr*nrp*rho*t2*t4*(5.0_r64*ic)-rn*nr*nrp*rho*t2*t6*(4.0_r64*ic) &
                +rn*nr*nrp*rhop*t2*t5*(8.0_r64*ic)+rn*nz*nzp*rho*t2*t5*ic-rn*nz*nzp*rhop*t2*t4*ic &
                -nr*nrp*rho*t2*t4*t10*(4.4e+1_r64*ic)+nr*nrp*rho*t2*t6*t10*(2.8e+1_r64*ic)-nr*nrp*rhop*t2*t5*t10*(3.2e+1_r64*ic) &
                -nz*nzp*rho*t2*t5*t10*(4.0_r64*ic)+nz*nzp*rhop*t2*t4*t10*(4.0_r64*ic)-rn*nr*nzp*t2*t5*zh*(4.0_r64*ic) &
                +rn*nrp*nz*t5*t12*zh*(4.0_r64*ic)+rn*nr*nzp*t4*t15*zh*ic+nr*nzp*t2*t5*t10*zh*(1.6e+1_r64*ic) &
                -nrp*nz*t5*t10*t12*zh*(1.6e+1_r64*ic)
                et27 = nr*nzp*t4*t10*t15*zh*-(4.0_r64*ic)-chi*rn*nr*nrp*rhop*t2*(2.0e+1_r64*ic)-chi*rn*nz*nzp*rho*t2*ic &
                +chi*rn*nz*nzp*rhop*t18*(4.0_r64*ic)+chi*nr*nrp*rhop*t2*t10*(8.0e+1_r64*ic)+chi*nz*nzp*rho*t2*t10*(4.0_r64*ic) &
                -chi*nz*nzp*rhop*t10*t18*(1.6e+1_r64*ic)+chi*rn*nr*nzp*t2*zh*(1.2e+1_r64*ic)-chi*rn*nrp*nz*t12*zh*(8.0_r64*ic) &
                -chi*rn*nrp*nz*t15*zh*(4.0_r64*ic)
                coE(4) = (t21*t23*t28*t31*t34**3*t50*(et25+et26+et27))/4.0_r64
                et28 = nr*nrp*t12*5.0_r64+nr*nrp*t15*5.0_r64+nz*nzp*t2*2.0_r64-chi*nr*nrp*t2*1.8e+1_r64+chi*nz*nzp*t18 &
                +nr*nrp*t2*t5*1.6e+1_r64-nr*nrp*t4*t12*4.0_r64-nr*nrp*t4*t15*4.0_r64+nr*nrp*t9*t12*1.6e+1_r64 &
                +nr*nrp*t9*t15*1.6e+1_r64-nz*nzp*t2*t4*2.0_r64+nz*nzp*t2*t27+nrp*nz*rho*zh*5.0_r64-nr*nzp*rhop*zh*5.0_r64 &
                -chi*nr*nrp*t2*t9*2.4e+1_r64-chi*nz*nzp*t9*t18*4.0_r64+chi*nr*nzp*rho*zh-chi*nrp*nz*rhop*zh &
                +nr*nrp*t2*t5*t9*3.2e+1_r64-nr*nrp*t4*t9*t12*2.0e+1_r64-nr*nrp*t4*t9*t15*2.0e+1_r64-nz*nzp*t2*t4*t9*4.0_r64 &
                -nrp*nz*rho*t4*zh*4.0_r64+nrp*nz*rho*t9*zh*1.6e+1_r64+nr*nzp*rhop*t4*zh*4.0_r64-nr*nzp*rhop*t9*zh*1.6e+1_r64 &
                -chi*nr*nzp*rho*t9*zh*4.0_r64
                et29 = chi*nrp*nz*rhop*t27*zh-nrp*nz*rho*t4*t9*zh*2.0e+1_r64+nr*nzp*rhop*t4*t9*zh*2.0e+1_r64
                coK(5) = (t21*t34**3*t38*(et28+et29))/(rho*8.0_r64)
                et30 = nr*nrp*t2*6.0_r64+nz*nzp*t18*3.0_r64+chi*nr*nrp*t12*8.0_r64+chi*nr*nrp*t15*8.0_r64+chi*nz*nzp*t2*2.0_r64 &
                -nr*nrp*t2*t4*3.0e+1_r64+nr*nrp*t2*t6*1.6e+1_r64-nr*nrp*t2*t9*1.6e+1_r64+nr*nrp*t2*t11*1.6e+1_r64 &
                -nr*nrp*t5*t12*4.0_r64-nr*nrp*t5*t15*4.0_r64-nz*nzp*t2*t5*2.0_r64+nz*nzp*t4*t18-nz*nzp*t9*t18*8.0_r64 &
                -nz*nzp*t11*t18*1.6e+1_r64+nr*nzp*rho*zh*3.0_r64-nrp*nz*rhop*zh*3.0_r64-chi*nr*nrp*t9*t12*2.8e+1_r64 &
                -chi*nr*nrp*t11*t12*1.6e+1_r64-chi*nr*nrp*t9*t15*2.8e+1_r64-chi*nr*nrp*t11*t15*1.6e+1_r64 &
                -chi*nz*nzp*t2*t9*8.0_r64+chi*nrp*nz*rho*zh*8.0_r64-chi*nr*nzp*rhop*zh*8.0_r64+nr*nrp*t2*t4*t9*9.6e+1_r64 &
                -nr*nrp*t2*t6*t9*4.8e+1_r64
                et31 = nr*nrp*t2*t6*t11*-1.6e+1_r64+nr*nrp*t5*t9*t12*1.2e+1_r64+nr*nrp*t5*t11*t12*1.6e+1_r64 &
                +nr*nrp*t5*t9*t15*1.2e+1_r64+nr*nrp*t5*t11*t15*1.6e+1_r64+nz*nzp*t2*t5*t9*8.0_r64-nz*nzp*t4*t9*t18*8.0_r64 &
                +nz*nzp*t4*t11*t18*1.6e+1_r64+nr*nzp*rho*t4*zh-nrp*nz*rho*t5*zh*4.0_r64-nr*nzp*rho*t9*zh*8.0_r64 &
                -nr*nzp*rho*t11*zh*1.6e+1_r64-nrp*nz*rhop*t4*zh+nr*nzp*rhop*t5*zh*4.0_r64+nrp*nz*rhop*t9*zh*8.0_r64 &
                +nrp*nz*rhop*t11*zh*1.6e+1_r64-chi*nrp*nz*rho*t9*zh*2.8e+1_r64-chi*nrp*nz*rho*t11*zh*1.6e+1_r64 &
                +chi*nr*nzp*rhop*t9*zh*2.8e+1_r64+chi*nr*nzp*rhop*t11*zh*1.6e+1_r64-nr*nzp*rho*t4*t9*zh*8.0_r64 &
                +nrp*nz*rho*t5*t9*zh*1.2e+1_r64+nr*nzp*rho*t4*t11*zh*1.6e+1_r64
                et32 = nrp*nz*rho*t5*t11*zh*1.6e+1_r64+nrp*nz*rhop*t4*t9*zh*8.0_r64-nr*nzp*rhop*t5*t9*zh*1.2e+1_r64 &
                -nrp*nz*rhop*t4*t11*zh*1.6e+1_r64-nr*nzp*rhop*t5*t11*zh*1.6e+1_r64
                coE(5) = (rhop*t21*t28*t31*t35*t50*(et30+et31+et32))/8.0_r64
                coK(6) = rn*t21*t23*t34**3*t38*(t37+t42+t48+t49+t51+t52+t53+t60+t62+t66+t67+t69+t70+t74+t75+t76+t77+t78 &
                +nrp*nz*rho*t2*(3.0_r64*ic)-nrp*nz*rho*t2*t4*(3.0_r64*ic))*(-1.0_r64/4.0_r64)
                coE(6) = (rn*rhop*t21*t28*t31*t35*(t40+t41+t54+t58+t63+t64+t68+nrp*nz*rho*t18*(4.0_r64*ic)-nrp*nz*rhop*t2*ic &
                -nr*nzp*rhop*t18*(5.0_r64*ic)+nr*nrp*t12*zh*(4.0_r64*ic)+nr*nrp*t15*zh*(5.0_r64*ic)+chi*nrp*nz*rho*t2*ic &
                -nrp*nz*rho*t2*t5*ic+nrp*nz*rhop*t2*t4*ic+nr*nzp*rhop*t4*t18*ic-nr*nrp*t4*t15*zh*ic))/(rho*4.0_r64)
                et33 = t45-nz*nzp*rho*t19*5.0_r64+nrp*nz*t2*t18*1.0e+1_r64-nr*nzp*t12*t18*5.0_r64-nr*nzp*t15*t18*5.0_r64 &
                +nr*nrp*t16*zh*5.0_r64-nr*nzp*t4*1.0_r64/t34**4*6.0_r64+nr*nzp*t6*1.0_r64/t34**4*2.0_r64 &
                -nr*nzp*t9*1.0_r64/t34**4*4.0_r64+chi*nz*nzp*rhop*t19*8.0_r64-chi*nr*nzp*t2*t12-chi*nr*nzp*t2*t15 &
                +chi*nr*nzp*t2*t18*1.7e+1_r64-chi*nrp*nz*t12*t18*1.0e+1_r64-chi*nrp*nz*t15*t18*8.0_r64 &
                -chi*nr*nrp*t13*zh*1.0e+1_r64-nz*nzp*rho*t4*t19*3.0_r64+nz*nzp*rho*t9*t19*2.0e+1_r64+nr*nzp*t2*t5*t12 &
                +nr*nzp*t2*t5*t15+nrp*nz*t2*t4*t18*6.0_r64-nr*nzp*t2*t5*t18-nrp*nz*t2*t9*t18*4.0e+1_r64-nr*nzp*t4*t12*t18*3.0_r64 &
                +nrp*nz*t5*t12*t18*2.0_r64-nr*nzp*t4*t15*t18*3.0_r64
                et34 = nr*nzp*t9*t12*t18*2.0e+1_r64+nr*nzp*t9*t15*t18*2.0e+1_r64+nr*nrp*rho*t2*zh*1.9e+1_r64 &
                +nz*nzp*rhop*t2*zh*2.0_r64+nr*nrp*t5*t13*zh*2.0_r64+nr*nrp*t4*t16*zh*3.0_r64-nr*nrp*t9*t16*zh*2.0e+1_r64 &
                +nr*nzp*t6*t27*1.0_r64/t34**4-chi*nr*nrp*rhop*t2*zh*2.6e+1_r64-chi*nz*nzp*rho*t2*zh*2.0_r64 &
                +chi*nr*nrp*t9*t13*zh*4.0e+1_r64+nz*nzp*rho*t4*t9*t19*1.2e+1_r64-nr*nzp*t2*t5*t9*t12*4.0_r64 &
                -nr*nzp*t2*t5*t9*t15*4.0_r64-nrp*nz*t2*t4*t9*t18*2.4e+1_r64+nr*nzp*t4*t9*t12*t18*1.2e+1_r64 &
                -nrp*nz*t5*t9*t12*t18*8.0_r64+nr*nzp*t4*t9*t15*t18*1.2e+1_r64+nr*nzp*t2*t5*t18*t27+nr*nrp*rho*t2*t4*zh &
                +nr*nrp*rho*t2*t6*zh*4.0_r64-nr*nrp*rho*t2*t9*zh*4.0e+1_r64+nr*nrp*rhop*t2*t5*zh*2.0_r64
                et35 = nz*nzp*rho*t2*t5*zh*2.0_r64-nz*nzp*rhop*t2*t4*zh*2.0_r64-nz*nzp*rhop*t2*t9*zh*8.0_r64 &
                -nr*nrp*t5*t9*t13*zh*8.0_r64-nr*nrp*t4*t9*t16*zh*1.2e+1_r64-chi*nz*nzp*rhop*t9*t19*3.2e+1_r64 &
                -chi*nr*nzp*t2*t9*t18*6.8e+1_r64+chi*nrp*nz*t9*t12*t18*4.0e+1_r64+chi*nr*nzp*t2*t12*t27 &
                +chi*nrp*nz*t9*t15*t18*3.2e+1_r64+chi*nr*nzp*t2*t15*t27+chi*nr*nrp*rhop*t2*t9*zh*1.04e+2_r64 &
                +chi*nz*nzp*rho*t2*t9*zh*8.0_r64-nr*nrp*rho*t2*t4*t9*zh*7.6e+1_r64+nr*nrp*rho*t2*t6*t9*zh*2.0e+1_r64 &
                -nr*nrp*rhop*t2*t5*t9*zh*8.0_r64-nz*nzp*rho*t2*t5*t9*zh*8.0_r64+nz*nzp*rhop*t2*t4*t9*zh*8.0_r64
                coK(7) = rhop*t21*t36*t39*(et33+et34+et35)*(-1.0_r64/8.0_r64)
                et36 = t46+t6*t46+nz*nzp*rhop*t19*9.0_r64-nr*nzp*t2*t12*3.0_r64-nr*nzp*t2*t15*3.0_r64+nr*nzp*t2*t18*2.1e+1_r64 &
                -nrp*nz*t12*t18*1.5e+1_r64-nrp*nz*t15*t18*9.0_r64-nr*nrp*t13*zh*1.5e+1_r64+chi*nr*nzp*1.0_r64/t34**4*1.0e+1_r64 &
                -nrp*nz*t4*1.0_r64/t34**4*8.0_r64-nr*nzp*t5*1.0_r64/t34**4*1.2e+1_r64+nr*nzp*t7*1.0_r64/t34**4*2.0_r64 &
                -chi*nz*nzp*rho*t19*2.9e+1_r64+chi*nrp*nz*t2*t18*5.8e+1_r64-chi*nr*nzp*t12*t18*2.9e+1_r64 &
                -chi*nr*nzp*t15*t18*2.9e+1_r64+chi*nr*nrp*t16*zh*2.9e+1_r64-nz*nzp*rho*t5*t19*3.0_r64 &
                +nz*nzp*rhop*t4*t19*2.3e+1_r64-nz*nzp*rhop*t9*t19*4.0_r64+nr*nzp*t2*t4*t12*2.0_r64+nr*nzp*t2*t6*t12 &
                +nr*nzp*t2*t4*t15*2.0_r64+nr*nzp*t2*t6*t15+nr*nzp*t2*t4*t18*4.4e+1_r64
                et37 = nrp*nz*t2*t5*t18*6.0_r64-nr*nzp*t2*t6*t18-nr*nzp*t2*t9*t18*4.0_r64-nrp*nz*t4*t12*t18*1.9e+1_r64 &
                -nr*nzp*t5*t12*t18*3.0_r64+nrp*nz*t6*t12*t18*2.0_r64-nrp*nz*t4*t15*t18*2.3e+1_r64-nr*nzp*t5*t15*t18*3.0_r64 &
                +nrp*nz*t15*t18*t27-nr*nrp*rhop*t2*zh*3.3e+1_r64-nz*nzp*rho*t2*zh*6.0_r64-nr*nrp*t4*t13*zh*1.9e+1_r64 &
                +nr*nrp*t6*t13*zh*2.0_r64+nr*nrp*t5*t16*zh*3.0_r64+chi*nr*nrp*rho*t2*zh*9.1e+1_r64+chi*nz*nzp*rhop*t2*zh*8.0_r64 &
                -chi*nr*nrp*t9*t16*zh*4.0_r64-nz*nzp*rho*t5*t9*t19*4.0_r64+nz*nzp*rhop*t4*t19*t27+nrp*nz*t2*t5*t9*t18*8.0_r64 &
                -nr*nzp*t5*t9*t12*t18*4.0_r64-nrp*nz*t6*t9*t12*t18*4.0_r64-nrp*nz*t4*t9*t15*t18*4.0_r64 &
                -nr*nzp*t5*t9*t15*t18*4.0_r64
                et38 = nr*nzp*t2*t6*t18*t27+nrp*nz*t4*t12*t18*t27+nr*nrp*rho*t2*t5*zh+nr*nrp*rho*t2*t7*zh*4.0_r64 &
                -nr*nrp*rhop*t2*t4*zh*6.5e+1_r64+nr*nrp*rhop*t2*t6*zh*2.0_r64+nr*nrp*rhop*t2*t27*zh+nz*nzp*rho*t2*t4*zh*4.0_r64 &
                +nz*nzp*rho*t2*t6*zh*2.0_r64-nz*nzp*rhop*t2*t5*zh*8.0_r64-nr*nrp*t6*t9*t13*zh*4.0_r64+nr*nrp*t4*t13*t27*zh &
                +nr*nrp*t5*t16*t27*zh+chi*nz*nzp*rho*t19*t27-chi*nrp*nz*t2*t9*t18*8.0_r64+chi*nr*nzp*t12*t18*t27 &
                +chi*nr*nzp*t15*t18*t27-chi*nr*nrp*rho*t2*t9*zh*8.0_r64+nr*nrp*rho*t2*t5*t27*zh+nr*nrp*rho*t2*t7*t27*zh &
                -nr*nrp*rhop*t2*t6*t9*zh*8.0_r64+nr*nrp*rhop*t2*t4*t27*zh
                coE(7) = (rhop*t21*t30*t32*t36*(et36+et37+et38))/8.0_r64
                coK(8) = rn*t21*t35*t38*(t37+t42+t48+t49+t51+t52+t53+t60+t62+t66+t67+t69+t70+t74+t75+t76+t77+t78 &
                -nr*nzp*rhop*t2*(3.0_r64*ic)+nr*nzp*rhop*t2*t4*(3.0_r64*ic))*(-1.0_r64/4.0_r64)
                coE(8) = (rn*t21*t28*t31*t35*(t40+t41+t54+t58+t63+t64+t68+nr*nzp*rho*t2*ic+nrp*nz*rho*t18*(5.0_r64*ic) &
                -nr*nzp*rhop*t18*(4.0_r64*ic)+nr*nrp*t12*zh*(5.0_r64*ic)+nr*nrp*t15*zh*(4.0_r64*ic)-chi*nr*nzp*rhop*t2*ic &
                -nr*nzp*rho*t2*t4*ic-nrp*nz*rho*t4*t18*ic+nr*nzp*rhop*t2*t5*ic-nr*nrp*t4*t12*zh*ic))/4.0_r64
                et39 = nr*nrp*1.0_r64/t34**4*4.0_r64-chi*nz*nzp*t20*8.0_r64-nrp*nz*rho*t19*5.0_r64+nr*nzp*rhop*t19*5.0_r64 &
                -nr*nrp*t12*t18*5.0_r64-nr*nrp*t15*t18*5.0_r64-nz*nzp*t2*t18*4.0_r64-nr*nrp*t4*1.0_r64/t34**4*6.0_r64 &
                +nr*nrp*t6*1.0_r64/t34**4*2.0_r64-nr*nrp*t9*1.0_r64/t34**4*4.0_r64-chi*nr*nzp*rho*t19*8.0_r64 &
                +chi*nrp*nz*rhop*t19*8.0_r64-chi*nr*nrp*t2*t12-chi*nr*nrp*t2*t15+chi*nr*nrp*t2*t18*1.7e+1_r64 &
                +chi*nz*nzp*t9*t20*3.2e+1_r64-nrp*nz*rho*t4*t19*3.0_r64+nrp*nz*rho*t9*t19*2.0e+1_r64+nr*nzp*rhop*t4*t19*3.0_r64 &
                -nr*nzp*rhop*t9*t19*2.0e+1_r64+nr*nrp*t2*t5*t12+nr*nrp*t2*t5*t15-nr*nrp*t2*t5*t18-nr*nrp*t4*t12*t18*3.0_r64 &
                -nr*nrp*t4*t15*t18*3.0_r64
                et40 = nr*nrp*t9*t12*t18*2.0e+1_r64+nr*nrp*t9*t15*t18*2.0e+1_r64+nz*nzp*t2*t4*t18*4.0_r64 &
                +nz*nzp*t2*t9*t18*1.6e+1_r64-nr*nzp*rho*t2*zh*2.0_r64+nrp*nz*rhop*t2*zh*2.0_r64+nr*nrp*t6*t27*1.0_r64/t34**4 &
                -chi*nrp*nz*rho*t2*zh*2.0_r64+chi*nr*nzp*rhop*t2*zh*2.0_r64+nrp*nz*rho*t4*t9*t19*1.2e+1_r64 &
                -nr*nzp*rhop*t4*t9*t19*1.2e+1_r64-nr*nrp*t2*t5*t9*t12*4.0_r64-nr*nrp*t2*t5*t9*t15*4.0_r64 &
                +nr*nrp*t4*t9*t12*t18*1.2e+1_r64+nr*nrp*t4*t9*t15*t18*1.2e+1_r64+nr*nrp*t2*t5*t18*t27 &
                -nz*nzp*t2*t4*t9*t18*1.6e+1_r64+nr*nzp*rho*t2*t4*zh*2.0_r64+nrp*nz*rho*t2*t5*zh*2.0_r64 &
                +nr*nzp*rho*t2*t9*zh*8.0_r64-nrp*nz*rhop*t2*t4*zh*2.0_r64-nr*nzp*rhop*t2*t5*zh*2.0_r64 &
                -nrp*nz*rhop*t2*t9*zh*8.0_r64
                et41 = chi*nr*nzp*rho*t9*t19*3.2e+1_r64-chi*nrp*nz*rhop*t9*t19*3.2e+1_r64-chi*nr*nrp*t2*t9*t18*6.8e+1_r64 &
                +chi*nr*nrp*t2*t12*t27+chi*nr*nrp*t2*t15*t27+chi*nrp*nz*rho*t2*t9*zh*8.0_r64-chi*nr*nzp*rhop*t2*t9*zh*8.0_r64 &
                -nr*nzp*rho*t2*t4*t9*zh*8.0_r64-nrp*nz*rho*t2*t5*t9*zh*8.0_r64+nrp*nz*rhop*t2*t4*t9*zh*8.0_r64 &
                +nr*nzp*rhop*t2*t5*t9*zh*8.0_r64
                coK(9) = (t21*t35*t39*(et39+et40+et41)*(-1.0_r64/8.0_r64))/rho
                et42 = -t47+nz*nzp*t20*9.0_r64+nr*nzp*rho*t19*9.0_r64-nrp*nz*rhop*t19*9.0_r64+nr*nrp*t2*t12*3.0_r64 &
                +nr*nrp*t2*t15*3.0_r64-nr*nrp*t2*t18*2.1e+1_r64+nz*nzp*t4*t20*2.3e+1_r64-nz*nzp*t9*t20*4.0_r64 &
                -chi*nr*nrp*1.0_r64/t34**4*1.0e+1_r64+nr*nrp*t5*1.0_r64/t34**4*1.2e+1_r64-nr*nrp*t7*1.0_r64/t34**4*2.0_r64 &
                +nz*nzp*t4*1.0_r64/t34**4*8.0_r64-nz*nzp*t6*1.0_r64/t34**4*4.0_r64+chi*nrp*nz*rho*t19*2.9e+1_r64 &
                -chi*nr*nzp*rhop*t19*2.9e+1_r64+chi*nr*nrp*t12*t18*2.9e+1_r64+chi*nr*nrp*t15*t18*2.9e+1_r64 &
                +chi*nz*nzp*t2*t18*1.6e+1_r64+nr*nzp*rho*t4*t19*2.3e+1_r64+nrp*nz*rho*t5*t19*3.0_r64-nr*nzp*rho*t9*t19*4.0_r64 &
                -nrp*nz*rhop*t4*t19*2.3e+1_r64-nr*nzp*rhop*t5*t19*3.0_r64+nrp*nz*rhop*t19*t27
                et43 = nr*nrp*t2*t4*t12*-2.0_r64-nr*nrp*t2*t6*t12-nr*nrp*t2*t4*t15*2.0_r64-nr*nrp*t2*t6*t15 &
                -nr*nrp*t2*t4*t18*4.4e+1_r64+nr*nrp*t2*t6*t18+nr*nrp*t5*t12*t18*3.0_r64+nr*nrp*t5*t15*t18*3.0_r64 &
                +nr*nrp*t2*t18*t27-nz*nzp*t2*t5*t18*1.6e+1_r64+nz*nzp*t4*t20*t27+nrp*nz*rho*t2*zh*6.0_r64 &
                -nr*nzp*rhop*t2*zh*6.0_r64+chi*nr*nzp*rho*t2*zh*8.0_r64-chi*nrp*nz*rhop*t2*zh*8.0_r64+nr*nzp*rho*t4*t19*t27 &
                +nrp*nz*rho*t5*t19*t27-nrp*nz*rhop*t4*t9*t19*4.0_r64-nr*nzp*rhop*t5*t9*t19*4.0_r64-nr*nrp*t2*t6*t9*t18*4.0_r64 &
                +nr*nrp*t5*t12*t18*t27+nr*nrp*t5*t15*t18*t27-nrp*nz*rho*t2*t4*zh*4.0_r64-nr*nzp*rho*t2*t5*zh*8.0_r64 &
                -nrp*nz*rho*t2*t6*zh*2.0_r64
                et44 = nr*nzp*rhop*t2*t4*zh*4.0_r64+nrp*nz*rhop*t2*t5*zh*8.0_r64+nr*nzp*rhop*t2*t6*zh*2.0_r64 &
                -chi*nrp*nz*rho*t9*t19*4.0_r64+chi*nr*nzp*rhop*t19*t27-chi*nr*nrp*t9*t12*t18*4.0_r64 &
                -chi*nr*nrp*t9*t15*t18*4.0_r64
                coE(9) = rhop*t21*t30*t32*t36*(et42+et43+et44)*(-1.0_r64/8.0_r64)
              do e = 1, 9
                ff(jp,e) = (coK(e)*vk + coE(e)*ve)*wsp(jp)
              end do
              end block
            end do
            do e = 1, 9
              ar = (e-1)/3 + 1; b = mod(e-1,3) + 1
              do l = 1, p
                A((ar-1)*nt+nidx(i), (b-1)*nso+cols+l, md+1) = &
                  A((ar-1)*nt+nidx(i), (b-1)*nso+cols+l, md+1) + sum(ff(1:p,e)*Lc(:,l))
              end do
            end do
          end do
        end do
      end do
    end do
    deallocate(tin, joa, ci, zc, zcn, nidx, znear, znearn, ikq, ikc, C1, C2, C3, C4, C5, As, Ad, A1, A2, A3, A4, Gcq, Gpc, ff)
  end subroutine axissymstok_dlpn_blockmat_nmode_r64
end module axissymstok_specialquad_mod
