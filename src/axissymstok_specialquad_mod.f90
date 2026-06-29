module axissymstok_specialquad_mod
  use axistokes3d_mod, only: r64, gauss_r64, lagrange_interp_r64
  use axisym_modal_green_mod, only: modal_green_r64, modal_green_all_r64, modal_green_all_far_r64
  use specialquad_mod, only: sdspecialquad_r64
  use axissymstok_kernelsplit_mod, only: axissymstok_slp_coef_r64, axissymstok_slp_coef_nmode_r64, &
       axissymstok_slpn_coef_r64, axissymstok_slpn_coef_nmode_r64, &
       axissymstok_dlp_coef_r64, axissymstok_dlp_coef_nmode_r64, weight_setup_r64, &
       axissymstok_dlp_aziquad_r64, axissymstok_slpn_aziquad_r64, axissymstok_dlpn_coef_r64, &
       axissymstok_dlpn_coef_nmode_r64, axissymstok_dlpn_aziquad_r64, &
       axissymstok_slppres_coef_nmode_r64, axissymstok_dlppres_coef_nmode_r64
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
  public :: axissymstok_slppres_blockmat_nmode_r64
  public :: axissymstok_dlppres_blockmat_nmode_r64
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
          call modal_green_r64(chi, 0_8, vk, ve, Fn, An, dFn)
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
    ! Two-level with Ksub=1 (dyadic-pole inner mesh).  Outer COARSE panel pq: FAR targets -> naive on the COARSE
    ! nodes, straight into A.  NEAR targets -> per dyadic sub-panel of pq, inner-near = close-eval (targets near
    ! that sub-panel), inner-far = naive (targets far from it; well-separated in the dyadic mesh -> accurate),
    ! folded back via IPk.  A middle panel is ONE sub-panel == the coarse panel (no split, no cross-panel bug).
    integer(8),   intent(in)    :: nt, p, np, M, iside, iclosed
    complex(r64), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: sws(p*np), tpan(np+1), mu
    complex(r64), intent(inout) :: A(3*nt, 3*np*p, M+1)
    integer(8) :: q, pq, k, nso, i, j, jp, iq, ia, l, e, md, ar, b, ne, npa, nk, nkc, coff
    real(r64)  :: twopi, sumws, sumwsc, rho, rhop, zh, rr2, chi, ws, spd, rt, zti, vk, ve, rn, n2, fn2, ifn2
    real(r64)  :: muinv, ipi, rr, srt, rrt, cm1, cp1, icm1, t1, tN, tm, denom, rlo, rhi, tgi, split
    real(r64)  :: vka(0:M), vea(0:M), vkmat(p,0:M), vemat(p,0:M)
    real(r64)  :: sk1a(p),sk3a(p),se3a(p),sk5a(p),sk7a(p),se7a(p),sk9a(p),se9a(p),cse5a(p),prefa(p),p0a(p),p1a(p),wsa(p)
    complex(r64) :: c2a(p), c4a(p), c6a(p), c8a(p), cse2a(p)
    real(r64)  :: se1, se5
    complex(r64) :: sk2, se2, sk4, se4, sk6, sk8
    complex(r64) :: Be_rr(p), Be_rt(p), Be_rz(p), &
                    Be_tr(p), Be_tt(p), Be_tz(p), &
                    Be_zr(p), Be_zt(p), Be_zz(p)
    complex(r64) :: ic, za, zb, cc1, cc2, cc3
    real(r64)    :: tglp(p), wglp(p), Dp(p,p), tglq(2*p), wglq(2*p), Dq(2*p,2*p), IP2(2*p,p)
    real(r64)    :: tc(p,np), rk(p), re2(2), IPk(p,p), IPe2(2,p), IPqc(2*p,p), wsq(2*p), wsp(p)
    complex(r64) :: Yp(p), Ypb(p), Yq(2*p), dYq(2*p), dYp(p), nvq(2*p), wxpq(2*p), ec2(2), xc(p,np)
    real(r64),    allocatable :: tin(:)
    integer(8),   allocatable :: joa(:), nidx(:), ci(:)
    logical,      allocatable :: cl(:), ikc(:)
    complex(r64), allocatable :: znear(:), zc(:), C1(:,:,:), C2(:,:,:), C3(:,:,:), Gcq(:,:), Gpc(:,:)
    real(r64),    allocatable :: As(:,:), A1(:,:), A2(:,:), A3(:,:), A4(:,:)
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
    allocate(tin(np+1+256)); ne = np+1; tin(1:ne) = tpan(1:np+1)
    do while (tin(2) > t1)
      split = 0.5_r64*(tin(1)+tin(2)); do i = ne, 2, -1; tin(i+1) = tin(i); end do; tin(2) = split; ne = ne+1
    end do
    do while (tin(ne-1) < tN)
      split = 0.5_r64*(tin(ne)+tin(ne-1)); tin(ne+1) = tin(ne); tin(ne) = split; ne = ne+1
    end do
    npa = ne-1
    allocate(joa(npa))
    do k = 1, npa
      tm = 0.5_r64*(tin(k)+tin(k+1)); j = 1
      do i = 1, np
        if (tm >= tpan(i) .and. tm < tpan(i+1)) j = i
      end do
      joa(k) = j
    end do
    allocate(cl(nt), ikc(nt), nidx(nt), ci(nt), znear(nt), zc(nt))
    allocate(C1(3*nt,3*q,M+1), C2(3*nt,3*q,M+1), C3(3*nt,3*q,M+1))
    allocate(As(q,nt), Ad(q,nt), A1(q,nt), A2(q,nt), A3(q,nt), A4(q,nt), Gcq(nt,q), Gpc(nt,p))
    A = (0.0_r64,0.0_r64)
    do pq = 1, np                                                   ! ===== outer: COARSE panel pq =====
      coff = (pq-1)*p
      sumws = 0.0_r64
      do jp = 1, p
        sumws = sumws + sws(coff+jp)
      end do
      nk = 0
      do i = 1, nt
        cl(i) = (abs(tx(i)-sxlo(pq)) + abs(tx(i)-sxhi(pq))) < 1.5_r64*sumws
        if (cl(i)) then; nk = nk + 1; nidx(nk) = i; znear(nk) = tx(i); end if
      end do
      ! ---- FAR: naive on the COARSE nodes -> straight into A ----
      do jp = 1, p
        Yp(jp) = xc(jp,pq); wsp(jp) = sws(coff+jp)
      end do
      do i = 1, nt
        if (cl(i)) cycle
        rt = real(tx(i),r64); zti = aimag(tx(i))
        do jp = 1, p
          rhop = real(Yp(jp),r64); zh = zti - aimag(Yp(jp)); rho = rt
          rr2 = (rho-rhop)**2 + zh*zh; chi = 1.0_r64 + rr2/(2.0_r64*rho*rhop)
          call modal_green_all_far_r64(chi, M, vka, vea)
          vkmat(jp,0:M) = vka(0:M); vemat(jp,0:M) = vea(0:M)
          rr = rho*rhop; srt = 1.0_r64/sqrt(rr); rrt = srt*srt*srt
          cm1 = chi-1.0_r64; cp1 = chi+1.0_r64; icm1 = 1.0_r64/cm1; wsa(jp) = wsp(jp)*muinv
          sk1a(jp) = -(1.0_r64/8.0_r64)*rhop*ipi*rrt*(rho*rho+rhop*rhop-4.0_r64*chi*rr)
          sk3a(jp) =  (1.0_r64/8.0_r64)*rhop*rhop*ipi*rrt*zh
          se3a(jp) =  (1.0_r64/8.0_r64)*rhop*ipi*rrt*icm1*zh*(rho - chi*rhop)
          sk5a(jp) =  (1.0_r64/2.0_r64)*chi*rhop*srt*ipi
          sk7a(jp) = -(1.0_r64/8.0_r64)*srt*ipi*zh
          se7a(jp) = -(1.0_r64/8.0_r64)*rhop*ipi*rrt*icm1*zh*(rhop - chi*rho)
          sk9a(jp) =  (1.0_r64/4.0_r64)*rhop*srt*ipi
          se9a(jp) =  (1.0_r64/8.0_r64)*rhop*ipi*rrt*icm1*zh*zh
          c2a(jp)  = -(1.0_r64/4.0_r64)*srt*ipi*ic*(rho - chi*rhop)
          c4a(jp)  =  (1.0_r64/4.0_r64)*rhop*rhop*ipi*rrt*ic*(rhop - chi*rho)
          c6a(jp)  = -(1.0_r64/4.0_r64)*ic*rhop*rhop*ipi*rrt*zh
          c8a(jp)  = -(1.0_r64/4.0_r64)*ic*srt*ipi*zh
          cse2a(jp)= -(3.0_r64/4.0_r64)*rhop*srt*ipi*ic*cp1
          cse5a(jp)= -(1.0_r64/2.0_r64)*rhop*srt*ipi*cp1
          prefa(jp)= -(1.0_r64/8.0_r64)*rhop*ipi*rrt*icm1
          p0a(jp)  = 2.0_r64*rr + chi*rho*rho + chi*rhop*rhop - 4.0_r64*chi*chi*rr
          p1a(jp)  = 4.0_r64*rr + 4.0_r64*chi*chi*rr - 4.0_r64*chi*rho*rho - 4.0_r64*chi*rhop*rhop
        end do
        do md = 0, M
          rn = real(md,r64); n2 = rn*rn; fn2 = 4.0_r64*n2 - 1.0_r64; ifn2 = 1.0_r64/fn2
          do jp = 1, p
            vk = vkmat(jp,md); ve = vemat(jp,md); ws = wsa(jp)
            se1 = prefa(jp)*ifn2*(p0a(jp) + n2*p1a(jp))
            sk2 = c2a(jp)*rn;  se2 = cse2a(jp)*(rn*ifn2);  sk4 = c4a(jp)*rn;  se4 = -se2
            se5 = cse5a(jp)*((n2-1.0_r64)*ifn2);  sk6 = c6a(jp)*rn;  sk8 = c8a(jp)*rn
            Be_rr(jp) = (sk1a(jp)*vk + se1*ve)*ws;  Be_rt(jp) = (sk2*vk + se2*ve)*ws;  Be_rz(jp) = (sk3a(jp)*vk + se3a(jp)*ve)*ws
            Be_tr(jp) = (sk4*vk + se4*ve)*ws;       Be_tt(jp) = (sk5a(jp)*vk + se5*ve)*ws;  Be_tz(jp) = (sk6*vk)*ws
            Be_zr(jp) = (sk7a(jp)*vk + se7a(jp)*ve)*ws;  Be_zt(jp) = (sk8*vk)*ws;        Be_zz(jp) = (sk9a(jp)*vk + se9a(jp)*ve)*ws
          end do
          do l = 1, p
            A(0*nt+i, 0*nso+coff+l, md+1) = A(0*nt+i, 0*nso+coff+l, md+1) + Be_rr(l)
            A(0*nt+i, 1*nso+coff+l, md+1) = A(0*nt+i, 1*nso+coff+l, md+1) + Be_rt(l)
            A(0*nt+i, 2*nso+coff+l, md+1) = A(0*nt+i, 2*nso+coff+l, md+1) + Be_rz(l)
            A(1*nt+i, 0*nso+coff+l, md+1) = A(1*nt+i, 0*nso+coff+l, md+1) + Be_tr(l)
            A(1*nt+i, 1*nso+coff+l, md+1) = A(1*nt+i, 1*nso+coff+l, md+1) + Be_tt(l)
            A(1*nt+i, 2*nso+coff+l, md+1) = A(1*nt+i, 2*nso+coff+l, md+1) + Be_tz(l)
            A(2*nt+i, 0*nso+coff+l, md+1) = A(2*nt+i, 0*nso+coff+l, md+1) + Be_zr(l)
            A(2*nt+i, 1*nso+coff+l, md+1) = A(2*nt+i, 1*nso+coff+l, md+1) + Be_zt(l)
            A(2*nt+i, 2*nso+coff+l, md+1) = A(2*nt+i, 2*nso+coff+l, md+1) + Be_zz(l)
          end do
        end do
      end do
      if (nk == 0) cycle
      ! ---- NEAR: per dyadic sub-panel of pq, inner-near (close-eval) + inner-far (naive) ----
      do k = 1, npa
        if (joa(k) /= pq) cycle
        denom = tpan(pq+1) - tpan(pq)
        do i = 1, p
          tgi = tin(k) + (1.0_r64+tglp(i))/2.0_r64*(tin(k+1)-tin(k))
          rk(i) = (2.0_r64*tgi - (tpan(pq)+tpan(pq+1)))/denom
        end do
        call lagrange_interp_r64(p, tglp, p, rk, IPk)
        Ypb = matmul(IPk, xc(:,pq))
        rlo = (2.0_r64*tin(k)   - (tpan(pq)+tpan(pq+1)))/denom
        rhi = (2.0_r64*tin(k+1) - (tpan(pq)+tpan(pq+1)))/denom
        re2(1) = rlo; re2(2) = rhi
        call lagrange_interp_r64(p, tglp, 2_8, re2, IPe2)
        ec2 = matmul(IPe2, xc(:,pq)); za = ec2(1); zb = ec2(2)
        IPqc = matmul(IP2, IPk)
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
          ikc(i) = (abs(znear(i)-za) + abs(znear(i)-zb)) < 1.5_r64*sumwsc
          if (ikc(i)) then; nkc = nkc + 1; ci(nkc) = i; zc(nkc) = znear(i); end if
        end do
        ! inner-near: close-eval the nkc targets near sub-panel k, fold q -> coarse via IPqc
        if (nkc > 0) then
          call axissymstok_slp_coef_nmode_r64(nkc, zc(1:nkc), q, Yq, M, mu, &
               C1(1:3*nkc,1:3*q,:), C2(1:3*nkc,1:3*q,:), C3(1:3*nkc,1:3*q,:))
          call sdspecialquad_r64(nkc, zc(1:nkc), q, Yq, nvq, wxpq, za, zb, iside, &
               As(1:q,1:nkc), Ad(1:q,1:nkc), A1(1:q,1:nkc), A2(1:q,1:nkc), A3(1:q,1:nkc), A4(1:q,1:nkc))
          do md = 0, M
            do e = 1, 9
              ar = (e-1)/3 + 1; b = mod(e-1,3) + 1
              do iq = 1, q
                do ia = 1, nkc
                  cc1 = C1(3*(ia-1)+ar, 3*(iq-1)+b, md+1); cc2 = C2(3*(ia-1)+ar, 3*(iq-1)+b, md+1)
                  cc3 = C3(3*(ia-1)+ar, 3*(iq-1)+b, md+1)
                  Gcq(ia,iq) = twopi*cc1*As(iq,ia) + cc2*wsq(iq) &
                             + cmplx(twopi*real(Ad(iq,ia)*(cc3/nvq(iq)), r64), 0.0_r64, r64)
                end do
              end do
              Gpc(1:nkc,1:p) = matmul(Gcq(1:nkc,1:q), IPqc)
              select case (e)
              case (1); do ia=1,nkc; do l=1,p; A(0*nt+nidx(ci(ia)), 0*nso+coff+l, md+1) = A(0*nt+nidx(ci(ia)), 0*nso+coff+l, md+1) + Gpc(ia,l); end do; end do
              case (2); do ia=1,nkc; do l=1,p; A(0*nt+nidx(ci(ia)), 1*nso+coff+l, md+1) = A(0*nt+nidx(ci(ia)), 1*nso+coff+l, md+1) + Gpc(ia,l); end do; end do
              case (3); do ia=1,nkc; do l=1,p; A(0*nt+nidx(ci(ia)), 2*nso+coff+l, md+1) = A(0*nt+nidx(ci(ia)), 2*nso+coff+l, md+1) + Gpc(ia,l); end do; end do
              case (4); do ia=1,nkc; do l=1,p; A(1*nt+nidx(ci(ia)), 0*nso+coff+l, md+1) = A(1*nt+nidx(ci(ia)), 0*nso+coff+l, md+1) + Gpc(ia,l); end do; end do
              case (5); do ia=1,nkc; do l=1,p; A(1*nt+nidx(ci(ia)), 1*nso+coff+l, md+1) = A(1*nt+nidx(ci(ia)), 1*nso+coff+l, md+1) + Gpc(ia,l); end do; end do
              case (6); do ia=1,nkc; do l=1,p; A(1*nt+nidx(ci(ia)), 2*nso+coff+l, md+1) = A(1*nt+nidx(ci(ia)), 2*nso+coff+l, md+1) + Gpc(ia,l); end do; end do
              case (7); do ia=1,nkc; do l=1,p; A(2*nt+nidx(ci(ia)), 0*nso+coff+l, md+1) = A(2*nt+nidx(ci(ia)), 0*nso+coff+l, md+1) + Gpc(ia,l); end do; end do
              case (8); do ia=1,nkc; do l=1,p; A(2*nt+nidx(ci(ia)), 1*nso+coff+l, md+1) = A(2*nt+nidx(ci(ia)), 1*nso+coff+l, md+1) + Gpc(ia,l); end do; end do
              case (9); do ia=1,nkc; do l=1,p; A(2*nt+nidx(ci(ia)), 2*nso+coff+l, md+1) = A(2*nt+nidx(ci(ia)), 2*nso+coff+l, md+1) + Gpc(ia,l); end do; end do
              end select
            end do
          end do
        end if
        ! inner-far: naive on sub-panel k for the near targets that are far from k (dyadic -> well separated), fold via IPk
        do i = 1, nk
          if (ikc(i)) cycle
          rt = real(znear(i),r64); zti = aimag(znear(i))
          do jp = 1, p
            rhop = real(Ypb(jp),r64); zh = zti - aimag(Ypb(jp)); rho = rt
            rr2 = (rho-rhop)**2 + zh*zh; chi = 1.0_r64 + rr2/(2.0_r64*rho*rhop)
            call modal_green_all_far_r64(chi, M, vka, vea)
            vkmat(jp,0:M) = vka(0:M); vemat(jp,0:M) = vea(0:M)
          rr = rho*rhop; srt = 1.0_r64/sqrt(rr); rrt = srt*srt*srt
          cm1 = chi-1.0_r64; cp1 = chi+1.0_r64; icm1 = 1.0_r64/cm1; wsa(jp) = wsp(jp)*muinv
          sk1a(jp) = -(1.0_r64/8.0_r64)*rhop*ipi*rrt*(rho*rho+rhop*rhop-4.0_r64*chi*rr)
          sk3a(jp) =  (1.0_r64/8.0_r64)*rhop*rhop*ipi*rrt*zh
          se3a(jp) =  (1.0_r64/8.0_r64)*rhop*ipi*rrt*icm1*zh*(rho - chi*rhop)
          sk5a(jp) =  (1.0_r64/2.0_r64)*chi*rhop*srt*ipi
          sk7a(jp) = -(1.0_r64/8.0_r64)*srt*ipi*zh
          se7a(jp) = -(1.0_r64/8.0_r64)*rhop*ipi*rrt*icm1*zh*(rhop - chi*rho)
          sk9a(jp) =  (1.0_r64/4.0_r64)*rhop*srt*ipi
          se9a(jp) =  (1.0_r64/8.0_r64)*rhop*ipi*rrt*icm1*zh*zh
          c2a(jp)  = -(1.0_r64/4.0_r64)*srt*ipi*ic*(rho - chi*rhop)
          c4a(jp)  =  (1.0_r64/4.0_r64)*rhop*rhop*ipi*rrt*ic*(rhop - chi*rho)
          c6a(jp)  = -(1.0_r64/4.0_r64)*ic*rhop*rhop*ipi*rrt*zh
          c8a(jp)  = -(1.0_r64/4.0_r64)*ic*srt*ipi*zh
          cse2a(jp)= -(3.0_r64/4.0_r64)*rhop*srt*ipi*ic*cp1
          cse5a(jp)= -(1.0_r64/2.0_r64)*rhop*srt*ipi*cp1
          prefa(jp)= -(1.0_r64/8.0_r64)*rhop*ipi*rrt*icm1
          p0a(jp)  = 2.0_r64*rr + chi*rho*rho + chi*rhop*rhop - 4.0_r64*chi*chi*rr
          p1a(jp)  = 4.0_r64*rr + 4.0_r64*chi*chi*rr - 4.0_r64*chi*rho*rho - 4.0_r64*chi*rhop*rhop
          end do
          do md = 0, M
            rn = real(md,r64); n2 = rn*rn; fn2 = 4.0_r64*n2 - 1.0_r64; ifn2 = 1.0_r64/fn2
            do jp = 1, p
            vk = vkmat(jp,md); ve = vemat(jp,md); ws = wsa(jp)
            se1 = prefa(jp)*ifn2*(p0a(jp) + n2*p1a(jp))
            sk2 = c2a(jp)*rn;  se2 = cse2a(jp)*(rn*ifn2);  sk4 = c4a(jp)*rn;  se4 = -se2
            se5 = cse5a(jp)*((n2-1.0_r64)*ifn2);  sk6 = c6a(jp)*rn;  sk8 = c8a(jp)*rn
            Be_rr(jp) = (sk1a(jp)*vk + se1*ve)*ws;  Be_rt(jp) = (sk2*vk + se2*ve)*ws;  Be_rz(jp) = (sk3a(jp)*vk + se3a(jp)*ve)*ws
            Be_tr(jp) = (sk4*vk + se4*ve)*ws;       Be_tt(jp) = (sk5a(jp)*vk + se5*ve)*ws;  Be_tz(jp) = (sk6*vk)*ws
            Be_zr(jp) = (sk7a(jp)*vk + se7a(jp)*ve)*ws;  Be_zt(jp) = (sk8*vk)*ws;        Be_zz(jp) = (sk9a(jp)*vk + se9a(jp)*ve)*ws
            end do
            do l = 1, p
              A(0*nt+nidx(i), 0*nso+coff+l, md+1) = A(0*nt+nidx(i), 0*nso+coff+l, md+1) + sum(Be_rr(1:p)*IPk(:,l))
              A(0*nt+nidx(i), 1*nso+coff+l, md+1) = A(0*nt+nidx(i), 1*nso+coff+l, md+1) + sum(Be_rt(1:p)*IPk(:,l))
              A(0*nt+nidx(i), 2*nso+coff+l, md+1) = A(0*nt+nidx(i), 2*nso+coff+l, md+1) + sum(Be_rz(1:p)*IPk(:,l))
              A(1*nt+nidx(i), 0*nso+coff+l, md+1) = A(1*nt+nidx(i), 0*nso+coff+l, md+1) + sum(Be_tr(1:p)*IPk(:,l))
              A(1*nt+nidx(i), 1*nso+coff+l, md+1) = A(1*nt+nidx(i), 1*nso+coff+l, md+1) + sum(Be_tt(1:p)*IPk(:,l))
              A(1*nt+nidx(i), 2*nso+coff+l, md+1) = A(1*nt+nidx(i), 2*nso+coff+l, md+1) + sum(Be_tz(1:p)*IPk(:,l))
              A(2*nt+nidx(i), 0*nso+coff+l, md+1) = A(2*nt+nidx(i), 0*nso+coff+l, md+1) + sum(Be_zr(1:p)*IPk(:,l))
              A(2*nt+nidx(i), 1*nso+coff+l, md+1) = A(2*nt+nidx(i), 1*nso+coff+l, md+1) + sum(Be_zt(1:p)*IPk(:,l))
              A(2*nt+nidx(i), 2*nso+coff+l, md+1) = A(2*nt+nidx(i), 2*nso+coff+l, md+1) + sum(Be_zz(1:p)*IPk(:,l))
            end do
          end do
        end do
      end do
    end do
    deallocate(tin, joa, nidx, ci, znear, zc, cl, ikc, C1, C2, C3, As, Ad, A1, A2, A3, A4, Gcq, Gpc)
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
          call modal_green_r64(chi, 0_8, vk, ve, Fn, An, dFn)
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
    real(r64)  :: twopi, sumwso, sumwsc, rho, rhop, zh, rr2, chi, ws, spd, rt, zti, nr, nz, vk, ve, Fn, An, dFn, rn, n2, ifn2
    real(r64)  :: vka(0:M), vea(0:M), Fna(0:M), Ana(0:M), dFna(0:M), vkmat(p,0:M), vemat(p,0:M)
    complex(r64) :: coK(9), coE(9), pcs(27), pcsm(27,p)   ! pcs/pcsm = md-indep far coK/coE pieces (inlined per node, see the far blocks)
    complex(r64) :: Be_rr(p), Be_rt(p), Be_rz(p), &       ! 9 entries of the 3x3 traction block, per panel node
                    Be_tr(p), Be_tt(p), Be_tz(p), &
                    Be_zr(p), Be_zt(p), Be_zz(p)
    real(r64)  :: st1, stN, tm, denom, rlo, rhi, tgi, split
    complex(r64) :: ic, zac, zbc, cc1, cc2, cc3, cc4, Slog, Dval, Dz, Acau, Ahy
    real(r64)    :: tglp(p), wglp(p), Dp(p,p), tglq(2*p), wglq(2*p), Dq(2*p,2*p), IP2(2*p,p)
    real(r64)    :: tc(p,np), rk(p), re2(2), Lc(p,p), IPe2(2,p), IPqc(2*p,p), wsq(2*p), wsp(p)
    complex(r64) :: Ypb(p), Yq(2*p), dYq(2*p), dYp(p), nvq(2*p), wxpq(2*p), ec2(2), xc(p,np)
    real(r64),    allocatable :: tin(:)
    integer(8),   allocatable :: joa(:), ci(:), nidx(:)
    logical,      allocatable :: ikq(:), ikc(:)
    complex(r64), allocatable :: zc(:), zcn(:), znear(:), znearn(:), C1(:,:,:), C2(:,:,:), C3(:,:,:), C4(:,:,:), Gcq(:,:), Gpc(:,:)
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
    allocate(As(q,nt), Ad(q,nt), A1(q,nt), A2(q,nt), A3(q,nt), A4(q,nt), Gcq(nt,q), Gpc(nt,p))
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
          call modal_green_all_far_r64(chi, M, vka, vea)        ! all modes ONCE per node
          ! ---- md-INDEPENDENT pieces of the slpn far coK/coE coefficients (inlined, computed ONCE per node) ----
          ! coK linear in n2(=mode^2) with at most a 1/(4n2-1)=ifn2 factor; the do-md loop then assembles:
          !   coK diag(1,3,5,7,9)=p0+n2*p2   coK cross(2,4,6,8)=n*p1   coE(1,5)=ifn2*(p0+n2*p2)
          !   coE(2,4)=n*ifn2*(p1+n2*p3)     coE(3,7,9)=p0             coE(6,8)=n*p1
          ! pcs order: coK1_0 coK1_2 coE1_0 coE1_2 coK2_1 coE2_1 coE2_3 coK3_0 coK3_2 coE3_0
          !            coK4_1 coE4_1 coE4_3 coK5_0 coK5_2 coE5_0 coE5_2 coK6_1 coE6_1 coK7_0
          !            coK7_2 coE7_0 coK8_1 coE8_1 coK9_0 coK9_2 coE9_0   (gen /tmp/gen_slpn_pieces.m; 1/pi=ipi; uses outer ic)
          block
            real(r64) :: ipi, r15, r25
            real(r64) :: t2,t3,t4,t5,t6,t7,t8,t11,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t25,t26,t27,t28
            real(r64) :: t29,t30,t31,t32,t33,t34,t37,t39,t40,t41,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56
            complex(r64) :: t23,t24,t35,t36,t38,t42,t45,t46,t57
            ipi = 1.0_r64/acos(-1.0_r64)
            t2=chi*rho; t3=chi*rhop; t4=nr*rho; t5=rho*rhop; t6=nz*zh; t7=chi+1.0_r64; t8=chi*chi
            t11=rho*rho; t13=rhop*rhop; t14=t13*rhop; t15=zh*zh; t16=nr*rhop*3.0_r64; t18=chi-1.0_r64
            t17=nr*t3; t19=-t2; t20=-t3; t21=1.0_r64/t7; t22=t8-1.0_r64
            r15=1.0_r64/t5**1.5_r64; r25=1.0_r64/t5**2.5_r64
            t23=t4*ic; t24=t6*ic; t25=t4*t11*3.0_r64; t26=nr*t13*3.0_r64; t31=t6*t11*3.0_r64; t32=t6*t13*3.0_r64
            t33=1.0_r64/t18; t38=nr*rhop*t2*2.0_r64*ic; t39=t2*t2*t4; t41=rhop*t2*t6*10.0_r64
            t48=rhop*t2*t4*16.0_r64; t51=t2*t2*t6; t52=t3*t3*t6; t53=chi*t2*t3*t6*2.0_r64
            t27=-t17; t28=t13*t17*6.0_r64; t29=rho+t20; t30=rhop+t19; t34=t33*t33; t35=t17*ic
            t36=rho*t24; t40=t2*t17*2.0_r64; t42=rho*t23; t46=-t38; t47=t3*t3*t17*2.0_r64; t49=1.0_r64/t22
            t50=-t41; t54=-t48; t55=t2*t2*t17*4.0_r64; t37=-t28; t45=-t35; t56=t4+t6+t27; t57=t23+t24+t45
            pcs(1)  = ipi*rhop*r25*t49*(nr*t14*(-3.0_r64)-t4*t5*7.0_r64-t5*t6*4.0_r64+nr*rhop*t2*t2*4.0_r64 &
                      +nr*t2*t13*11.0_r64+rho*t2*t4+rho*t2*t6+rhop*t3*t6+rhop*t3*t17*2.0_r64+t2*t3*t6*2.0_r64 &
                      -t2*t3*t17*8.0_r64)*(-1.0_r64/8.0_r64)
            pcs(2)  = ipi*rhop*r25*t49*t56*(t5+rho*t19+rhop*t20+t2*t3)*(-1.0_r64/2.0_r64)
            pcs(3)  = ipi*rhop*t21*t34*r25*(t25+t31+t32+t37+t39+t47+t50+t51+t52+t53+t54+t55+t4*t13*3.0_r64 &
                      +rhop*t2*t17*17.0_r64-chi*t2*t3*t17*8.0_r64)*(-1.0_r64/8.0_r64)
            pcs(4)  = ipi*(rhop*t21*t34*r25*(t25+t31+t32+t37+t39+t47+t50+t51+t52+t53+t54+t55+t4*t13*6.0_r64 &
                      +rhop*t2*t17*11.0_r64-chi*t2*t3*t17*5.0_r64))/2.0_r64
            pcs(5)  = ipi*rhop*r15*(t4*2.0_r64*ic-t17*2.0_r64*ic+t24)*(-3.0_r64/4.0_r64)
            pcs(6)  = ipi*(t33*r15*(t36+t42+t46+nr*t13*3.0_r64*ic-t3*t6*ic-t3*t17*2.0_r64*ic))/4.0_r64
            pcs(7)  = ipi*(-t29*t33*r15*t57)
            pcs(8)  = ipi*rhop*r25*t49*zh*(t26+rho*t4+rho*t6-t3*t17*2.0_r64+t6*t20-nr*rhop*t2*2.0_r64)*(-1.0_r64/8.0_r64)
            pcs(9)  = ipi*(rhop*t29*r25*t49*t56*zh)/2.0_r64
            pcs(10) = ipi*rhop*t21*t34*r25*zh*(t40+rhop*t4*6.0_r64+rhop*t6*3.0_r64-rhop*t17*6.0_r64-t2*t4*4.0_r64 &
                      -t2*t6*4.0_r64+chi*t3*t6+chi*t3*t17*2.0_r64)*(-1.0_r64/8.0_r64)
            pcs(11) = ipi*t13*r25*(t36+t42+t46+nr*t13*ic)*(3.0_r64/4.0_r64)
            pcs(12) = ipi*(t13*t33*r25*(rhop*t4*(-4.0_r64*ic)-rhop*t6*ic+rhop*t35+t2*t17*2.0_r64*ic+t2*t23+t2*t24))/4.0_r64
            pcs(13) = ipi*t13*t30*t33*r25*t57
            pcs(14) = ipi*rhop*r15*(t4+t6-t17*4.0_r64)*(-1.0_r64/4.0_r64)
            pcs(15) = ipi*rhop*r15*t56*(-1.0_r64/2.0_r64)
            pcs(16) = ipi*rhop*t33*r15*(t16+chi*t6-chi*t17*4.0_r64+nr*t2)*(-1.0_r64/4.0_r64)
            pcs(17) = ipi*(rhop*t33*r15*(t16+chi*t6*2.0_r64-chi*t17*5.0_r64+nr*t2*2.0_r64))/2.0_r64
            pcs(18) = ipi*nr*t14*r25*zh*cmplx(0.0_r64,-0.75_r64,r64)
            pcs(19) = ipi*t13*t33*r25*t57*zh*(-1.0_r64/4.0_r64)
            pcs(20) = ipi*rhop*r25*t49*zh*(t40-rhop*t4*4.0_r64-rhop*t6+rhop*t17+t2*t4+t2*t6)*(-1.0_r64/8.0_r64)
            pcs(21) = ipi*rhop*t30*r25*t49*t56*zh*(-1.0_r64/2.0_r64)
            pcs(22) = ipi*(rhop*t21*t34*r25*zh*(t26+chi*t40+rho*t4*3.0_r64+rho*t6*3.0_r64-t3*t6*4.0_r64+t3*t17 &
                      +nr*t2*t2+chi*t2*t6-nr*rhop*t2*10.0_r64))/8.0_r64
            pcs(23) = ipi*t4*t13*r25*zh*cmplx(0.0_r64,-0.75_r64,r64)
            pcs(24) = ipi*t33*r15*t57*zh*(-1.0_r64/4.0_r64)
            pcs(25) = ipi*rhop*t15*r25*t49*t56*(-1.0_r64/8.0_r64)
            pcs(26) = ipi*(rhop*t15*r25*t49*t56)/2.0_r64
            pcs(27) = ipi*rhop*t15*t21*t34*r25*(t16-chi*t6*4.0_r64+chi*t17-nr*t2*4.0_r64)*(-1.0_r64/8.0_r64)
          end block
          do md = 0, M
            rn = real(md,r64); n2 = rn*rn; ifn2 = 1.0_r64/(4.0_r64*n2-1.0_r64)
            vk = vka(md); ve = vea(md)
            coK(1)=pcs(1)+n2*pcs(2);  coK(3)=pcs(8)+n2*pcs(9);  coK(5)=pcs(14)+n2*pcs(15)
            coK(7)=pcs(20)+n2*pcs(21); coK(9)=pcs(25)+n2*pcs(26)
            coK(2)=rn*pcs(5); coK(4)=rn*pcs(11); coK(6)=rn*pcs(18); coK(8)=rn*pcs(23)
            coE(1)=ifn2*(pcs(3)+n2*pcs(4)); coE(5)=ifn2*(pcs(16)+n2*pcs(17))
            coE(2)=rn*ifn2*(pcs(6)+n2*pcs(7)); coE(4)=rn*ifn2*(pcs(12)+n2*pcs(13))
            coE(3)=pcs(10); coE(7)=pcs(22); coE(9)=pcs(27)
            coE(6)=rn*pcs(19); coE(8)=rn*pcs(24)
            ! 3x3 traction block entries  Be_** = -(coK*vk + coE*ve)*ws   (row r/t/z = 0/1/2*nt, col r/t/z = 0/1/2*nso)
            Be_rr(jp) = -(coK(1)*vk+coE(1)*ve)*ws;  Be_rt(jp) = -(coK(2)*vk+coE(2)*ve)*ws;  Be_rz(jp) = -(coK(3)*vk+coE(3)*ve)*ws
            Be_tr(jp) = -(coK(4)*vk+coE(4)*ve)*ws;  Be_tt(jp) = -(coK(5)*vk+coE(5)*ve)*ws;  Be_tz(jp) = -(coK(6)*vk+coE(6)*ve)*ws
            Be_zr(jp) = -(coK(7)*vk+coE(7)*ve)*ws;  Be_zt(jp) = -(coK(8)*vk+coE(8)*ve)*ws;  Be_zz(jp) = -(coK(9)*vk+coE(9)*ve)*ws
            A(0*nt+i, 0*nso+jj, md+1) = Be_rr(jp);  A(0*nt+i, 1*nso+jj, md+1) = Be_rt(jp);  A(0*nt+i, 2*nso+jj, md+1) = Be_rz(jp)
            A(1*nt+i, 0*nso+jj, md+1) = Be_tr(jp);  A(1*nt+i, 1*nso+jj, md+1) = Be_tt(jp);  A(1*nt+i, 2*nso+jj, md+1) = Be_tz(jp)
            A(2*nt+i, 0*nso+jj, md+1) = Be_zr(jp);  A(2*nt+i, 1*nso+jj, md+1) = Be_zt(jp);  A(2*nt+i, 2*nso+jj, md+1) = Be_zz(jp)
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
              select case (e)                                  ! fold Gpc into the (row,col) 3x3 block of A
              case (1); do ia=1,nkc; do l=1,p; A(0*nt+nidx(ci(ia)), 0*nso+cols+l, md+1) = A(0*nt+nidx(ci(ia)), 0*nso+cols+l, md+1) + Gpc(ia,l); end do; end do  ! Be_rr
              case (2); do ia=1,nkc; do l=1,p; A(0*nt+nidx(ci(ia)), 1*nso+cols+l, md+1) = A(0*nt+nidx(ci(ia)), 1*nso+cols+l, md+1) + Gpc(ia,l); end do; end do  ! Be_rt
              case (3); do ia=1,nkc; do l=1,p; A(0*nt+nidx(ci(ia)), 2*nso+cols+l, md+1) = A(0*nt+nidx(ci(ia)), 2*nso+cols+l, md+1) + Gpc(ia,l); end do; end do  ! Be_rz
              case (4); do ia=1,nkc; do l=1,p; A(1*nt+nidx(ci(ia)), 0*nso+cols+l, md+1) = A(1*nt+nidx(ci(ia)), 0*nso+cols+l, md+1) + Gpc(ia,l); end do; end do  ! Be_tr
              case (5); do ia=1,nkc; do l=1,p; A(1*nt+nidx(ci(ia)), 1*nso+cols+l, md+1) = A(1*nt+nidx(ci(ia)), 1*nso+cols+l, md+1) + Gpc(ia,l); end do; end do  ! Be_tt
              case (6); do ia=1,nkc; do l=1,p; A(1*nt+nidx(ci(ia)), 2*nso+cols+l, md+1) = A(1*nt+nidx(ci(ia)), 2*nso+cols+l, md+1) + Gpc(ia,l); end do; end do  ! Be_tz
              case (7); do ia=1,nkc; do l=1,p; A(2*nt+nidx(ci(ia)), 0*nso+cols+l, md+1) = A(2*nt+nidx(ci(ia)), 0*nso+cols+l, md+1) + Gpc(ia,l); end do; end do  ! Be_zr
              case (8); do ia=1,nkc; do l=1,p; A(2*nt+nidx(ci(ia)), 1*nso+cols+l, md+1) = A(2*nt+nidx(ci(ia)), 1*nso+cols+l, md+1) + Gpc(ia,l); end do; end do  ! Be_zt
              case (9); do ia=1,nkc; do l=1,p; A(2*nt+nidx(ci(ia)), 2*nso+cols+l, md+1) = A(2*nt+nidx(ci(ia)), 2*nso+cols+l, md+1) + Gpc(ia,l); end do; end do  ! Be_zz
              end select
            end do
          end do
        end if
        ! ---- inner far: naive slpn9far (negated) on the FINER subpanel nodes, TARGET normals, fold via Lc ----
        do i = 1, nk
          if (ikc(i)) cycle
          rt = real(znear(i),r64); zti = aimag(znear(i)); nr = real(znearn(i),r64); nz = aimag(znearn(i))
          do jp = 1, p                                          ! carrier_all ONCE per node (was M+1 modal_green_r64 calls)
            rhop = real(Ypb(jp),r64); zh = zti - aimag(Ypb(jp)); rho = rt
            rr2 = (rho-rhop)**2 + zh*zh; chi = 1.0_r64 + rr2/(2.0_r64*rho*rhop)
            call modal_green_all_far_r64(chi, M, vka, vea)
            vkmat(jp,0:M) = vka(0:M); vemat(jp,0:M) = vea(0:M)
            ! ---- md-indep far coK/coE pieces, inlined ONCE per node (see outer-far block for the pcs layout) ----
            block
              real(r64) :: ipi, r15, r25
              real(r64) :: t2,t3,t4,t5,t6,t7,t8,t11,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t25,t26,t27,t28
              real(r64) :: t29,t30,t31,t32,t33,t34,t37,t39,t40,t41,t47,t48,t49,t50,t51,t52,t53,t54,t55,t56
              complex(r64) :: t23,t24,t35,t36,t38,t42,t45,t46,t57
              ipi = 1.0_r64/acos(-1.0_r64)
              t2=chi*rho; t3=chi*rhop; t4=nr*rho; t5=rho*rhop; t6=nz*zh; t7=chi+1.0_r64; t8=chi*chi
              t11=rho*rho; t13=rhop*rhop; t14=t13*rhop; t15=zh*zh; t16=nr*rhop*3.0_r64; t18=chi-1.0_r64
              t17=nr*t3; t19=-t2; t20=-t3; t21=1.0_r64/t7; t22=t8-1.0_r64
              r15=1.0_r64/t5**1.5_r64; r25=1.0_r64/t5**2.5_r64
              t23=t4*ic; t24=t6*ic; t25=t4*t11*3.0_r64; t26=nr*t13*3.0_r64; t31=t6*t11*3.0_r64; t32=t6*t13*3.0_r64
              t33=1.0_r64/t18; t38=nr*rhop*t2*2.0_r64*ic; t39=t2*t2*t4; t41=rhop*t2*t6*10.0_r64
              t48=rhop*t2*t4*16.0_r64; t51=t2*t2*t6; t52=t3*t3*t6; t53=chi*t2*t3*t6*2.0_r64
              t27=-t17; t28=t13*t17*6.0_r64; t29=rho+t20; t30=rhop+t19; t34=t33*t33; t35=t17*ic
              t36=rho*t24; t40=t2*t17*2.0_r64; t42=rho*t23; t46=-t38; t47=t3*t3*t17*2.0_r64; t49=1.0_r64/t22
              t50=-t41; t54=-t48; t55=t2*t2*t17*4.0_r64; t37=-t28; t45=-t35; t56=t4+t6+t27; t57=t23+t24+t45
              pcsm(1,jp)  = ipi*rhop*r25*t49*(nr*t14*(-3.0_r64)-t4*t5*7.0_r64-t5*t6*4.0_r64+nr*rhop*t2*t2*4.0_r64 &
                        +nr*t2*t13*11.0_r64+rho*t2*t4+rho*t2*t6+rhop*t3*t6+rhop*t3*t17*2.0_r64+t2*t3*t6*2.0_r64 &
                        -t2*t3*t17*8.0_r64)*(-1.0_r64/8.0_r64)
              pcsm(2,jp)  = ipi*rhop*r25*t49*t56*(t5+rho*t19+rhop*t20+t2*t3)*(-1.0_r64/2.0_r64)
              pcsm(3,jp)  = ipi*rhop*t21*t34*r25*(t25+t31+t32+t37+t39+t47+t50+t51+t52+t53+t54+t55+t4*t13*3.0_r64 &
                        +rhop*t2*t17*17.0_r64-chi*t2*t3*t17*8.0_r64)*(-1.0_r64/8.0_r64)
              pcsm(4,jp)  = ipi*(rhop*t21*t34*r25*(t25+t31+t32+t37+t39+t47+t50+t51+t52+t53+t54+t55+t4*t13*6.0_r64 &
                        +rhop*t2*t17*11.0_r64-chi*t2*t3*t17*5.0_r64))/2.0_r64
              pcsm(5,jp)  = ipi*rhop*r15*(t4*2.0_r64*ic-t17*2.0_r64*ic+t24)*(-3.0_r64/4.0_r64)
              pcsm(6,jp)  = ipi*(t33*r15*(t36+t42+t46+nr*t13*3.0_r64*ic-t3*t6*ic-t3*t17*2.0_r64*ic))/4.0_r64
              pcsm(7,jp)  = ipi*(-t29*t33*r15*t57)
              pcsm(8,jp)  = ipi*rhop*r25*t49*zh*(t26+rho*t4+rho*t6-t3*t17*2.0_r64+t6*t20-nr*rhop*t2*2.0_r64)*(-1.0_r64/8.0_r64)
              pcsm(9,jp)  = ipi*(rhop*t29*r25*t49*t56*zh)/2.0_r64
              pcsm(10,jp) = ipi*rhop*t21*t34*r25*zh*(t40+rhop*t4*6.0_r64+rhop*t6*3.0_r64-rhop*t17*6.0_r64-t2*t4*4.0_r64 &
                        -t2*t6*4.0_r64+chi*t3*t6+chi*t3*t17*2.0_r64)*(-1.0_r64/8.0_r64)
              pcsm(11,jp) = ipi*t13*r25*(t36+t42+t46+nr*t13*ic)*(3.0_r64/4.0_r64)
              pcsm(12,jp) = ipi*(t13*t33*r25*(rhop*t4*(-4.0_r64*ic)-rhop*t6*ic+rhop*t35+t2*t17*2.0_r64*ic+t2*t23+t2*t24))/4.0_r64
              pcsm(13,jp) = ipi*t13*t30*t33*r25*t57
              pcsm(14,jp) = ipi*rhop*r15*(t4+t6-t17*4.0_r64)*(-1.0_r64/4.0_r64)
              pcsm(15,jp) = ipi*rhop*r15*t56*(-1.0_r64/2.0_r64)
              pcsm(16,jp) = ipi*rhop*t33*r15*(t16+chi*t6-chi*t17*4.0_r64+nr*t2)*(-1.0_r64/4.0_r64)
              pcsm(17,jp) = ipi*(rhop*t33*r15*(t16+chi*t6*2.0_r64-chi*t17*5.0_r64+nr*t2*2.0_r64))/2.0_r64
              pcsm(18,jp) = ipi*nr*t14*r25*zh*cmplx(0.0_r64,-0.75_r64,r64)
              pcsm(19,jp) = ipi*t13*t33*r25*t57*zh*(-1.0_r64/4.0_r64)
              pcsm(20,jp) = ipi*rhop*r25*t49*zh*(t40-rhop*t4*4.0_r64-rhop*t6+rhop*t17+t2*t4+t2*t6)*(-1.0_r64/8.0_r64)
              pcsm(21,jp) = ipi*rhop*t30*r25*t49*t56*zh*(-1.0_r64/2.0_r64)
              pcsm(22,jp) = ipi*(rhop*t21*t34*r25*zh*(t26+chi*t40+rho*t4*3.0_r64+rho*t6*3.0_r64-t3*t6*4.0_r64+t3*t17 &
                        +nr*t2*t2+chi*t2*t6-nr*rhop*t2*10.0_r64))/8.0_r64
              pcsm(23,jp) = ipi*t4*t13*r25*zh*cmplx(0.0_r64,-0.75_r64,r64)
              pcsm(24,jp) = ipi*t33*r15*t57*zh*(-1.0_r64/4.0_r64)
              pcsm(25,jp) = ipi*rhop*t15*r25*t49*t56*(-1.0_r64/8.0_r64)
              pcsm(26,jp) = ipi*(rhop*t15*r25*t49*t56)/2.0_r64
              pcsm(27,jp) = ipi*rhop*t15*t21*t34*r25*(t16-chi*t6*4.0_r64+chi*t17-nr*t2*4.0_r64)*(-1.0_r64/8.0_r64)
            end block
          end do
          do md = 0, M
            rn = real(md,r64); n2 = rn*rn; ifn2 = 1.0_r64/(4.0_r64*n2-1.0_r64)
            do jp = 1, p
              vk = vkmat(jp,md); ve = vemat(jp,md)
              coK(1)=pcsm(1,jp)+n2*pcsm(2,jp);   coK(3)=pcsm(8,jp)+n2*pcsm(9,jp);   coK(5)=pcsm(14,jp)+n2*pcsm(15,jp)
              coK(7)=pcsm(20,jp)+n2*pcsm(21,jp); coK(9)=pcsm(25,jp)+n2*pcsm(26,jp)
              coK(2)=rn*pcsm(5,jp); coK(4)=rn*pcsm(11,jp); coK(6)=rn*pcsm(18,jp); coK(8)=rn*pcsm(23,jp)
              coE(1)=ifn2*(pcsm(3,jp)+n2*pcsm(4,jp)); coE(5)=ifn2*(pcsm(16,jp)+n2*pcsm(17,jp))
              coE(2)=rn*ifn2*(pcsm(6,jp)+n2*pcsm(7,jp)); coE(4)=rn*ifn2*(pcsm(12,jp)+n2*pcsm(13,jp))
              coE(3)=pcsm(10,jp); coE(7)=pcsm(22,jp); coE(9)=pcsm(27,jp)
              coE(6)=rn*pcsm(19,jp); coE(8)=rn*pcsm(24,jp)
              ! 3x3 traction block entries  Be_** = -(coK*vk + coE*ve)*wsp   (per subpanel node jp, folded below)
              Be_rr(jp) = -(coK(1)*vk+coE(1)*ve)*wsp(jp);  Be_rt(jp) = -(coK(2)*vk+coE(2)*ve)*wsp(jp);  Be_rz(jp) = -(coK(3)*vk+coE(3)*ve)*wsp(jp)
              Be_tr(jp) = -(coK(4)*vk+coE(4)*ve)*wsp(jp);  Be_tt(jp) = -(coK(5)*vk+coE(5)*ve)*wsp(jp);  Be_tz(jp) = -(coK(6)*vk+coE(6)*ve)*wsp(jp)
              Be_zr(jp) = -(coK(7)*vk+coE(7)*ve)*wsp(jp);  Be_zt(jp) = -(coK(8)*vk+coE(8)*ve)*wsp(jp);  Be_zz(jp) = -(coK(9)*vk+coE(9)*ve)*wsp(jp)
            end do
            do l = 1, p                                        ! fold each named block (sum over subpanel nodes via Lc) into A
              A(0*nt+nidx(i), 0*nso+cols+l, md+1) = A(0*nt+nidx(i), 0*nso+cols+l, md+1) + sum(Be_rr(1:p)*Lc(:,l))   ! Be_rr
              A(0*nt+nidx(i), 1*nso+cols+l, md+1) = A(0*nt+nidx(i), 1*nso+cols+l, md+1) + sum(Be_rt(1:p)*Lc(:,l))   ! Be_rt
              A(0*nt+nidx(i), 2*nso+cols+l, md+1) = A(0*nt+nidx(i), 2*nso+cols+l, md+1) + sum(Be_rz(1:p)*Lc(:,l))   ! Be_rz
              A(1*nt+nidx(i), 0*nso+cols+l, md+1) = A(1*nt+nidx(i), 0*nso+cols+l, md+1) + sum(Be_tr(1:p)*Lc(:,l))   ! Be_tr
              A(1*nt+nidx(i), 1*nso+cols+l, md+1) = A(1*nt+nidx(i), 1*nso+cols+l, md+1) + sum(Be_tt(1:p)*Lc(:,l))   ! Be_tt
              A(1*nt+nidx(i), 2*nso+cols+l, md+1) = A(1*nt+nidx(i), 2*nso+cols+l, md+1) + sum(Be_tz(1:p)*Lc(:,l))   ! Be_tz
              A(2*nt+nidx(i), 0*nso+cols+l, md+1) = A(2*nt+nidx(i), 0*nso+cols+l, md+1) + sum(Be_zr(1:p)*Lc(:,l))   ! Be_zr
              A(2*nt+nidx(i), 1*nso+cols+l, md+1) = A(2*nt+nidx(i), 1*nso+cols+l, md+1) + sum(Be_zt(1:p)*Lc(:,l))   ! Be_zt
              A(2*nt+nidx(i), 2*nso+cols+l, md+1) = A(2*nt+nidx(i), 2*nso+cols+l, md+1) + sum(Be_zz(1:p)*Lc(:,l))   ! Be_zz
            end do
          end do
        end do
      end do
    end do
    deallocate(tin, joa, ci, zc, zcn, nidx, znear, znearn, ikq, ikc, C1, C2, C3, C4, As, Ad, A1, A2, A3, A4, Gcq, Gpc)
  end subroutine axissymstok_slpn_blockmat_nmode_r64

  ! ==================================================================================================
  ! Mode-n axisymmetric Stokes SINGLE-LAYER PRESSURE block (mu-free), per source force [f_rho,f_phi,f_z].
  !   A(nt, 3*np*p, M+1) : pressure p(x) = A * [f_rho; f_phi; f_z] of the SLP field.  Two-level close-eval
  !   (mirror axissymstok_slpn_blockmat_nmode_r64): outer 1.85 -> Ksub(+dyadic) subpanels -> inner 1.5.
  !   Near = the (private) axissymstok_slppres_coef_nmode_r64 split (Log+Smooth+Cauchy) on the 2p upsample;
  !   far  = the inlined carrier coefficients PK*VK + PE*VE (the leading carrier multipliers, like slpn far).
  ! ==================================================================================================
  subroutine axissymstok_slppres_blockmat_nmode_r64(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, A)
    integer(8),   intent(in)    :: nt, p, np, M, iside, iclosed
    complex(r64), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: sws(p*np), tpan(np+1)
    complex(r64), intent(inout) :: A(nt, 3*np*p, M+1)
    integer(8), parameter :: Ksub = 4
    integer(8) :: q, nso, i, j, jp, iq, ia, l, b, c, pq, qo, ne, npa, nk, nkc, jj, cols, md
    real(r64)  :: twopi, sumwso, sumwsc, rho, rhop, zh, rr2, chi, ws, spd, rt, zti, vk, ve, Fn, An, dFn
    real(r64)  :: vka(0:M), vea(0:M), Fna(0:M), Ana(0:M), dFna(0:M), vkmat(p,0:M), vemat(p,0:M)
    real(r64)  :: st1, stN, tm, denom, rlo, rhi, tgi, split
    complex(r64) :: ic, zac, zbc, Slog, Dval, cc1, cc2, cc3
    real(r64)    :: tglp(p), wglp(p), Dp(p,p), tglq(2*p), wglq(2*p), Dq(2*p,2*p), IP2(2*p,p)
    real(r64)    :: tc(p,np), rk(p), re2(2), Lc(p,p), IPe2(2,p), IPqc(2*p,p), wsq(2*p), wsp(p)
    complex(r64) :: Ypb(p), Yq(2*p), dYq(2*p), dYp(p), nvq(2*p), wxpq(2*p), ec2(2), xc(p,np)
    real(r64),    allocatable :: tin(:)
    integer(8),   allocatable :: joa(:), ci(:), nidx(:)
    logical,      allocatable :: ikq(:), ikc(:)
    complex(r64), allocatable :: zc(:), znear(:), Gcq(:,:), Gpc(:,:), ff(:,:), Ad(:,:), C1(:,:,:), C2(:,:,:), C3(:,:,:)
    real(r64),    allocatable :: As(:,:), A1(:,:), A2(:,:), A3(:,:), A4(:,:)
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
    allocate(As(q,nt), Ad(q,nt), A1(q,nt), A2(q,nt), A3(q,nt), A4(q,nt), Gcq(nt,q), Gpc(nt,p), ff(p,3))
    allocate(C1(nt,3*q,M+1), C2(nt,3*q,M+1), C3(nt,3*q,M+1))
    A = (0.0_r64,0.0_r64)
    do pq = 1, np                                                   ! ===== outer level: original panel pq =====
      cols = (pq-1)*p
      sumwso = 0.0_r64
      do jp = 1, p; sumwso = sumwso + sws(cols+jp); end do
      nk = 0
      do i = 1, nt
        ikq(i) = (abs(tx(i)-sxlo(pq)) + abs(tx(i)-sxhi(pq))) < 1.85_r64*sumwso
        if (ikq(i)) then; nk = nk + 1; nidx(nk) = i; znear(nk) = tx(i); end if
      end do
      ! ---- outer far: inlined carrier coefficients PK*VK + PE*VE on the coarse original nodes ----
      do i = 1, nt
        if (ikq(i)) cycle
        rt = real(tx(i),r64); zti = aimag(tx(i)); rho = rt
        do jp = 1, p
          jj = cols + jp; rhop = real(sx(jj),r64); zh = zti - aimag(sx(jj)); ws = sws(jj)
          rr2 = (rho-rhop)**2 + zh*zh; chi = 1.0_r64 + rr2/(2.0_r64*rho*rhop)
          call modal_green_all_r64(chi, M, vka, vea, Fna, Ana, dFna)        ! all modes ONCE per node
          do md = 0, M
            vk = vka(md); ve = vea(md)
            block
              real(r64) :: t2, t4, t6, t7, t8, rnf
              complex(r64) :: PKf(3), PEf(3)
              t2=rho*rhop; t4=1.0_r64/acos(-1.0_r64); t6=1.0_r64/(chi-1.0_r64)
              t7=1.0_r64/t2**1.5_r64; t8=(t4/sqrt(t2))/4.0_r64; rnf=real(md,r64)
              PKf(1)=cmplx(-t8,0.0_r64,r64); PEf(1)=cmplx(rhop*t4*t6*t7*(rhop-chi*rho)*(-0.25_r64),0.0_r64,r64)
              PKf(2)=rnf*t4/sqrt(t2)*(-0.5_r64)*ic; PEf(2)=(0.0_r64,0.0_r64)
              PKf(3)=(0.0_r64,0.0_r64); PEf(3)=cmplx(rhop*t4*t6*t7*zh/4.0_r64,0.0_r64,r64)
              do b = 1, 3
                A(i, (b-1)*nso+jj, md+1) = (PKf(b)*vk + PEf(b)*ve)*ws
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
        do jp = 1, p; wsp(jp) = wglp(jp)*abs(dYp(jp)); end do
        sumwsc = sum(wsp)
        nkc = 0
        do i = 1, nk
          ikc(i) = (abs(znear(i)-zac) + abs(znear(i)-zbc)) < 1.5_r64*sumwsc
          if (ikc(i)) then; nkc = nkc + 1; ci(nkc) = i; zc(nkc) = znear(i); end if
        end do
        ! ---- inner near: (private) pressure split coef -> Slog/Dval assembly, fold q->p_sub (IPqc) ----
        if (nkc > 0) then
          call axissymstok_slppres_coef_nmode_r64(nkc, zc(1:nkc), q, Yq, M, &
               C1(1:nkc,1:3*q,:), C2(1:nkc,1:3*q,:), C3(1:nkc,1:3*q,:))
          call sdspecialquad_r64(nkc, zc(1:nkc), q, Yq, nvq, wxpq, zac, zbc, iside, &
               As(1:q,1:nkc), Ad(1:q,1:nkc), A1(1:q,1:nkc), A2(1:q,1:nkc), A3(1:q,1:nkc), A4(1:q,1:nkc))
          do md = 0, M
            do b = 1, 3
              do iq = 1, q
                do ia = 1, nkc
                  cc1 = C1(ia, 3*(iq-1)+b, md+1); cc2 = C2(ia, 3*(iq-1)+b, md+1); cc3 = C3(ia, 3*(iq-1)+b, md+1)
                  Slog = cmplx(As(iq,ia),0.0_r64,r64); Dval = Ad(iq,ia)
                  Gcq(ia,iq) = twopi*cc1*Slog + cc2*wsq(iq) + cmplx(twopi*real(Dval*(cc3/nvq(iq)),r64),0.0_r64,r64)
                end do
              end do
              Gpc(1:nkc,1:p) = matmul(Gcq(1:nkc,1:q), IPqc)
              do ia = 1, nkc
                do l = 1, p
                  A(nidx(ci(ia)), (b-1)*nso+cols+l, md+1) = A(nidx(ci(ia)), (b-1)*nso+cols+l, md+1) + Gpc(ia,l)
                end do
              end do
            end do
          end do
        end if
        ! ---- inner far: inlined carrier coefficients on the FINER subpanel nodes, fold via Lc ----
        do i = 1, nk
          if (ikc(i)) cycle
          rt = real(znear(i),r64); zti = aimag(znear(i)); rho = rt
          do jp = 1, p                                          ! carrier_all ONCE per node (was M+1 modal_green_r64 calls)
            rhop = real(Ypb(jp),r64); zh = zti - aimag(Ypb(jp))
            rr2 = (rho-rhop)**2 + zh*zh; chi = 1.0_r64 + rr2/(2.0_r64*rho*rhop)
            call modal_green_all_r64(chi, M, vka, vea, Fna, Ana, dFna)
            vkmat(jp,0:M) = vka(0:M); vemat(jp,0:M) = vea(0:M)
          end do
          do md = 0, M
            do jp = 1, p
              rhop = real(Ypb(jp),r64); zh = zti - aimag(Ypb(jp))
              rr2 = (rho-rhop)**2 + zh*zh; chi = 1.0_r64 + rr2/(2.0_r64*rho*rhop)
              vk = vkmat(jp,md); ve = vemat(jp,md)
              block
                real(r64) :: t2, t4, t6, t7, t8, rnf
                complex(r64) :: PKf(3), PEf(3)
                t2=rho*rhop; t4=1.0_r64/acos(-1.0_r64); t6=1.0_r64/(chi-1.0_r64)
                t7=1.0_r64/t2**1.5_r64; t8=(t4/sqrt(t2))/4.0_r64; rnf=real(md,r64)
                PKf(1)=cmplx(-t8,0.0_r64,r64); PEf(1)=cmplx(rhop*t4*t6*t7*(rhop-chi*rho)*(-0.25_r64),0.0_r64,r64)
                PKf(2)=rnf*t4/sqrt(t2)*(-0.5_r64)*ic; PEf(2)=(0.0_r64,0.0_r64)
                PKf(3)=(0.0_r64,0.0_r64); PEf(3)=cmplx(rhop*t4*t6*t7*zh/4.0_r64,0.0_r64,r64)
                do b = 1, 3
                  ff(jp,b) = (PKf(b)*vk + PEf(b)*ve)*wsp(jp)
                end do
              end block
            end do
            do b = 1, 3
              do l = 1, p
                A(nidx(i), (b-1)*nso+cols+l, md+1) = A(nidx(i), (b-1)*nso+cols+l, md+1) + sum(ff(1:p,b)*Lc(:,l))
              end do
            end do
          end do
        end do
      end do
    end do
    deallocate(tin, joa, ci, zc, nidx, znear, ikq, ikc, As, Ad, A1, A2, A3, A4, Gcq, Gpc, ff, C1, C2, C3)
  end subroutine axissymstok_slppres_blockmat_nmode_r64

  ! ==================================================================================================
  ! Mode-n axisymmetric Stokes DOUBLE-LAYER (stresslet) PRESSURE block (mu-free), per source density.
  !   A(nt, 3*np*p, M+1) : pressure p^D(x) = A * [mu_rho; mu_phi; mu_z].  Two-level close-eval (mirror dlpn):
  !   near = the (private) axissymstok_dlppres_coef_nmode_r64 split (Log+Smooth+Cauchy+deriv-Cauchy) on the
  !   2p upsample; far = a carrier-free DIRECT azimuthal ring integral (nv pts) -- avoids carrier overflow.
  ! ==================================================================================================
  subroutine axissymstok_dlppres_blockmat_nmode_r64(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, A)
    integer(8),   intent(in)    :: nt, p, np, M, iside, iclosed
    complex(r64), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: sws(p*np), tpan(np+1)
    complex(r64), intent(inout) :: A(nt, 3*np*p, M+1)
    integer(8), parameter :: Ksub = 4, nv = 160
    integer(8) :: q, nso, i, j, jp, iq, ia, l, b, c, pq, qo, ne, npa, nk, nkc, jj, cols, md, iv
    real(r64)  :: twopi, sumwso, sumwsc, rho, rhop, zh, rr2, chi, ws, spd, rt, zti
    real(r64)  :: st1, stN, tm, denom, rlo, rhi, tgi, split, nrs, nzs, cv, sv, r2v, ndr, rcv, fac, vphi
    complex(r64) :: ic, zac, zbc, Slog, Dval, Dz, cc1, cc2, cc3, cc4, cq, hq, ev, ir3, ir5
    real(r64)    :: tglp(p), wglp(p), Dp(p,p), tglq(2*p), wglq(2*p), Dq(2*p,2*p), IP2(2*p,p)
    real(r64)    :: tc(p,np), rk(p), re2(2), Lc(p,p), IPe2(2,p), IPqc(2*p,p), wsq(2*p), wsp(p)
    complex(r64) :: Ypb(p), Yq(2*p), dYq(2*p), dYp(p), nvq(2*p), nYp(p), wxpq(2*p), ec2(2), xc(p,np)
    complex(r64) :: Pr(M+1), Pp(M+1), Pz(M+1)
    real(r64),    allocatable :: tin(:)
    integer(8),   allocatable :: joa(:), ci(:), nidx(:)
    logical,      allocatable :: ikq(:), ikc(:)
    complex(r64), allocatable :: zc(:), znear(:), Gcq(:,:), Gpc(:,:), Ad(:,:), C1(:,:,:), C2(:,:,:), C3(:,:,:), C4(:,:,:)
    real(r64),    allocatable :: As(:,:), A1(:,:), A2(:,:), A3(:,:), A4(:,:)
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
    allocate(As(q,nt), Ad(q,nt), A1(q,nt), A2(q,nt), A3(q,nt), A4(q,nt), Gcq(nt,q), Gpc(nt,p))
    allocate(C1(nt,3*q,M+1), C2(nt,3*q,M+1), C3(nt,3*q,M+1), C4(nt,3*q,M+1))
    A = (0.0_r64,0.0_r64)
    do pq = 1, np                                                   ! ===== outer level: original panel pq =====
      cols = (pq-1)*p
      sumwso = 0.0_r64
      do jp = 1, p; sumwso = sumwso + sws(cols+jp); end do
      nk = 0
      do i = 1, nt
        ikq(i) = (abs(tx(i)-sxlo(pq)) + abs(tx(i)-sxhi(pq))) < 1.85_r64*sumwso
        if (ikq(i)) then; nk = nk + 1; nidx(nk) = i; znear(nk) = tx(i); end if
      end do
      ! ---- outer far: carrier-free azimuthal ring integral on the coarse original nodes ----
      do i = 1, nt
        if (ikq(i)) cycle
        rt = real(tx(i),r64); zti = aimag(tx(i)); rho = rt
        do jp = 1, p
          jj = cols + jp; rhop = real(sx(jj),r64); zh = zti - aimag(sx(jj)); ws = sws(jj)
          nrs = real(snx(jj),r64); nzs = aimag(snx(jj))
          Pr = (0.0_r64,0.0_r64); Pp = (0.0_r64,0.0_r64); Pz = (0.0_r64,0.0_r64)
          do iv = 1, nv
            vphi = twopi*real(iv-1,r64)/real(nv,r64); cv = cos(vphi); sv = sin(vphi)
            r2v = rho*rho + rhop*rhop - 2.0_r64*rho*rhop*cv + zh*zh
            rcv = rho*cv - rhop; ndr = nrs*rcv + nzs*zh
            do md = 0, M
              ev = exp(ic*real(md,r64)*vphi); ir3 = ev/r2v**1.5_r64; ir5 = ev/r2v**2.5_r64
              Pr(md+1) = Pr(md+1) + (-nrs*ir3 + 3.0_r64*ndr*rcv*ir5)
              Pp(md+1) = Pp(md+1) + (3.0_r64*ndr*(-rho*sv)*ir5)
              Pz(md+1) = Pz(md+1) + (-nzs*ir3 + 3.0_r64*ndr*zh*ir5)
            end do
          end do
          fac = rhop*ws/real(nv,r64)
          do md = 0, M
            A(i,         jj, md+1) = Pr(md+1)*fac
            A(i,   nso+jj, md+1)   = Pp(md+1)*fac
            A(i, 2*nso+jj, md+1)   = Pz(md+1)*fac
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
        do jp = 1, p; wsp(jp) = wglp(jp)*abs(dYp(jp)); nYp(jp) = -ic*dYp(jp)/abs(dYp(jp)); end do
        sumwsc = sum(wsp)
        nkc = 0
        do i = 1, nk
          ikc(i) = (abs(znear(i)-zac) + abs(znear(i)-zbc)) < 1.5_r64*sumwsc
          if (ikc(i)) then; nkc = nkc + 1; ci(nkc) = i; zc(nkc) = znear(i); end if
        end do
        ! ---- inner near: (private) stresslet-pressure split coef -> 4-bucket assembly, fold q->p_sub ----
        if (nkc > 0) then
          call axissymstok_dlppres_coef_nmode_r64(nkc, zc(1:nkc), q, Yq, nvq, M, &
               C1(1:nkc,1:3*q,:), C2(1:nkc,1:3*q,:), C3(1:nkc,1:3*q,:), C4(1:nkc,1:3*q,:))
          call sdspecialquad_r64(nkc, zc(1:nkc), q, Yq, nvq, wxpq, zac, zbc, iside, &
               As(1:q,1:nkc), Ad(1:q,1:nkc), A1(1:q,1:nkc), A2(1:q,1:nkc), A3(1:q,1:nkc), A4(1:q,1:nkc))
          do md = 0, M
            do b = 1, 3
              do iq = 1, q
                do ia = 1, nkc
                  cc1 = C1(ia,3*(iq-1)+b,md+1); cc2 = C2(ia,3*(iq-1)+b,md+1)
                  cc3 = C3(ia,3*(iq-1)+b,md+1); cc4 = C4(ia,3*(iq-1)+b,md+1)
                  Slog = cmplx(As(iq,ia),0.0_r64,r64); Dval = Ad(iq,ia); Dz = cmplx(A1(iq,ia),-A2(iq,ia),r64)
                  if (b == 2) then                                  ! phi-force entry: imaginary partner (Re -> i Im)
                    cq =  ic*twopi*aimag(Dval*(cc3/nvq(iq)))
                    hq = -ic*twopi*aimag(Dz  *(cc4/nvq(iq)))
                  else
                    cq = cmplx( twopi*real(Dval*(cc3/nvq(iq)), r64), 0.0_r64, r64)
                    hq = cmplx(-twopi*real(Dz  *(cc4/nvq(iq)), r64), 0.0_r64, r64)
                  end if
                  Gcq(ia,iq) = twopi*cc1*Slog + cc2*wsq(iq) + cq + hq
                end do
              end do
              Gpc(1:nkc,1:p) = matmul(Gcq(1:nkc,1:q), IPqc)
              do ia = 1, nkc
                do l = 1, p
                  A(nidx(ci(ia)), (b-1)*nso+cols+l, md+1) = A(nidx(ci(ia)), (b-1)*nso+cols+l, md+1) + Gpc(ia,l)
                end do
              end do
            end do
          end do
        end if
        ! ---- inner far: ring integral on the FINER subpanel nodes, fold via Lc ----
        do i = 1, nk
          if (ikc(i)) cycle
          rt = real(znear(i),r64); zti = aimag(znear(i)); rho = rt
          do jp = 1, p
            rhop = real(Ypb(jp),r64); zh = zti - aimag(Ypb(jp)); ws = wsp(jp)
            nrs = real(nYp(jp),r64); nzs = aimag(nYp(jp))
            Pr = (0.0_r64,0.0_r64); Pp = (0.0_r64,0.0_r64); Pz = (0.0_r64,0.0_r64)
            do iv = 1, nv
              vphi = twopi*real(iv-1,r64)/real(nv,r64); cv = cos(vphi); sv = sin(vphi)
              r2v = rho*rho + rhop*rhop - 2.0_r64*rho*rhop*cv + zh*zh
              rcv = rho*cv - rhop; ndr = nrs*rcv + nzs*zh
              do md = 0, M
                ev = exp(ic*real(md,r64)*vphi); ir3 = ev/r2v**1.5_r64; ir5 = ev/r2v**2.5_r64
                Pr(md+1) = Pr(md+1) + (-nrs*ir3 + 3.0_r64*ndr*rcv*ir5)
                Pp(md+1) = Pp(md+1) + (3.0_r64*ndr*(-rho*sv)*ir5)
                Pz(md+1) = Pz(md+1) + (-nzs*ir3 + 3.0_r64*ndr*zh*ir5)
              end do
            end do
            fac = rhop*ws/real(nv,r64)
            do md = 0, M
              do l = 1, p
                A(nidx(i),         cols+l, md+1) = A(nidx(i),         cols+l, md+1) + Pr(md+1)*fac*Lc(jp,l)
                A(nidx(i),   nso+cols+l, md+1)   = A(nidx(i),   nso+cols+l, md+1)   + Pp(md+1)*fac*Lc(jp,l)
                A(nidx(i), 2*nso+cols+l, md+1)   = A(nidx(i), 2*nso+cols+l, md+1)   + Pz(md+1)*fac*Lc(jp,l)
              end do
            end do
          end do
        end do
      end do
    end do
    deallocate(tin, joa, ci, zc, nidx, znear, ikq, ikc, As, Ad, A1, A2, A3, A4, Gcq, Gpc, C1, C2, C3, C4)
  end subroutine axissymstok_dlppres_blockmat_nmode_r64

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
          call modal_green_r64(chi, 0_8, vk, ve, Fn, An, dFn)
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
    real(r64)  :: vka(0:M), vea(0:M), Fna(0:M), Ana(0:M), dFna(0:M), vkmat(p,0:M), vemat(p,0:M)
    complex(r64) :: coK(9), coE(9), pcs(27), pcsm(27,p)
    complex(r64) :: Be_rr(p), Be_rt(p), Be_rz(p), &       ! 9 entries of the 3x3 block, per panel node
                    Be_tr(p), Be_tt(p), Be_tz(p), &
                    Be_zr(p), Be_zt(p), Be_zz(p)
    real(r64)  :: n2, ifn2
    real(r64)  :: st1, stN, tm, denom, rlo, rhi, tgi, split
    complex(r64) :: ic, zac, zbc, cc1, cc2, cc3, cc4, Slog, Dval, Dz, Acau, Ahy
    real(r64)    :: tglp(p), wglp(p), Dp(p,p), tglq(2*p), wglq(2*p), Dq(2*p,2*p), IP2(2*p,p)
    real(r64)    :: tc(p,np), rk(p), re2(2), Lc(p,p), IPe2(2,p), IPqc(2*p,p), wsq(2*p), wsp(p)
    complex(r64) :: Ypc(p), Ypb(p), Yq(2*p), dYq(2*p), dYp(p), nvq(2*p), wxpq(2*p), ec2(2), xc(p,np)
    real(r64),    allocatable :: tin(:)
    integer(8),   allocatable :: joa(:), ci(:), nidx(:)
    logical,      allocatable :: ikq(:), ikc(:)
    complex(r64), allocatable :: zc(:), znear(:), C1(:,:,:), C2(:,:,:), C3(:,:,:), C4(:,:,:), Gcq(:,:), Gpc(:,:)
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
    allocate(As(q,nt), Ad(q,nt), A1(q,nt), A2(q,nt), A3(q,nt), A4(q,nt), Gcq(nt,q), Gpc(nt,p))
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
          call modal_green_all_far_r64(chi, M, vka, vea)        ! all modes ONCE per node
          block
            real(r64) :: ipi, r15, r25
            real(r64) :: t2,t3,t4,t5,t6,t7,t8,t11,t12,t13,t15,t16,t17,t18,t19,t20,t21,t22,t23,t26,t27,t28,t29,t30,t31,t32
            real(r64) :: t33,t34,t36,t39,t40,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55
            complex(r64) :: t24,t25,t35,t37,t41,t56
            ipi = 1.0_r64/acos(-1.0_r64)
            t2=chi*rho; t3=chi*rhop; t4=nrp*rhop; t5=rho*rhop; t6=nzp*zh; t7=chi+1.0_r64; t8=chi*chi
            t11=rho*rho; t12=t11*rho; t13=rhop*rhop; t15=zh*zh; t16=nrp*rho*3.0_r64; t18=chi-1.0_r64
            t17=nrp*t2; t19=-t2; t20=-t3; t21=-t4; t22=1.0_r64/t7; t23=t8-1.0_r64
            r15=1.0_r64/t5**1.5_r64; r25=1.0_r64/t5**2.5_r64
            t24=t4*ic; t25=t6*ic; t26=nrp*t11*3.0_r64; t27=t4*t13*3.0_r64; t31=t6*t11*3.0_r64; t32=t6*t13*3.0_r64
            t33=1.0_r64/t18; t40=rhop*t2*t6*10.0_r64; t41=nrp*t11*3.0_r64*ic; t45=rhop*t2*t4*16.0_r64
            t48=t2*t2*t6; t49=t3*t3*t6; t50=chi*t2*t3*t6*2.0_r64
            t28=t11*t17*6.0_r64; t29=rho+t20; t30=rhop+t19; t34=t33*t33; t35=-t24; t36=-t27; t37=t17*ic
            t39=t3*t17*2.0_r64; t44=t2*t2*t17*2.0_r64; t46=1.0_r64/t23; t47=-t40; t52=t3*t4*t20; t53=t3*t3*t17*4.0_r64
            t55=t6+t17+t21; t51=-t44; t54=-t53; t56=t25+t35+t37
            pcs(1) = ipi*rhop*r25*t46*(nrp*t12*3.0_r64+t4*t5*7.0_r64-t5*t6*4.0_r64-rho*t2*t4*11.0_r64+rho*t2*t6-rho*t2*t17*2.0_r64 &
                    +rhop*t3*t6+rhop*t4*t20-t2*t3*t4*4.0_r64+t2*t3*t6*2.0_r64+t2*t3*t17*8.0_r64)*(-1.0_r64/8.0_r64)
            pcs(2) = ipi*rhop*r25*t46*t55*(t5+rho*t19+rhop*t20+t2*t3)*(-1.0_r64/2.0_r64)
            pcs(3) = ipi*rhop*t22*t34*r25*(t28+t31+t32+t36+t45+t47+t48+t49+t50+t51+t52+t54-t4*t11*3.0_r64 &
                    -t2*t2*t4*17.0_r64+chi*t2*t3*t17*8.0_r64)*(-1.0_r64/8.0_r64)
            pcs(4) = ipi*(rhop*t22*t34*r25*(t28+t31+t32+t36+t45+t47+t48+t49+t50+t51+t52+t54-t4*t11*6.0_r64 &
                    -t2*t2*t4*11.0_r64+chi*t2*t3*t17*5.0_r64))/2.0_r64
            pcs(5) = ipi*(r15*(t41+rhop*t4*3.0_r64*ic-rhop*t6*3.0_r64*ic-t2*t4*6.0_r64*ic))/4.0_r64
            pcs(6) = ipi*(t33*r15*(rho*t4*(-4.0_r64*ic)+rho*t25+rho*t37-t3*t6*ic+t3*t17*2.0_r64*ic+t3*t24))/4.0_r64
            pcs(7) = ipi*(-t29*t33*r15*t56)
            pcs(8) = ipi*rhop*r25*t46*zh*(t39-rho*t4*4.0_r64+rho*t6+rho*t17+t3*t4+t6*t20)*(-1.0_r64/8.0_r64)
            pcs(9) = ipi*(rhop*t29*r25*t46*t55*zh)/2.0_r64
            pcs(10) = ipi*(rhop*t22*t34*r25*zh*(t26+chi*t39+rhop*t4*3.0_r64-rhop*t6*3.0_r64-t2*t4*10.0_r64+t2*t6*4.0_r64 &
                    +t2*t17+nrp*t3*t3+chi*t6*t20))/8.0_r64
            pcs(11) = ipi*rhop*r15*(t4*(-2.0_r64*ic)+t17*2.0_r64*ic+t25)*(3.0_r64/4.0_r64)
            pcs(12) = ipi*(t13*t33*r25*(t41-rhop*t6*ic+rhop*t24-t2*t4*2.0_r64*ic-t2*t17*2.0_r64*ic+t2*t25))/4.0_r64
            pcs(13) = ipi*t13*t30*t33*r25*t56
            pcs(14) = ipi*rhop*r15*(t17*3.0_r64+t55)*(-1.0_r64/4.0_r64)
            pcs(15) = ipi*rhop*r15*t55*(-1.0_r64/2.0_r64)
            pcs(16) = ipi*(rhop*t33*r15*(t16-chi*t6-chi*t17*4.0_r64+nrp*t3))/4.0_r64
            pcs(17) = ipi*rhop*t33*r15*(t16-chi*t6*2.0_r64-chi*t17*5.0_r64+nrp*t3*2.0_r64)*(-1.0_r64/2.0_r64)
            pcs(18) = ipi*t4*r15*zh*cmplx(0.0_r64,0.75_r64,r64)
            pcs(19) = ipi*t13*t33*r25*t56*zh*(-1.0_r64/4.0_r64)
            pcs(20) = ipi*rhop*r25*t46*zh*(t26+rhop*t4-rhop*t6-t2*t4*2.0_r64+t2*t6-t2*t17*2.0_r64)*(-1.0_r64/8.0_r64)
            pcs(21) = ipi*rhop*t30*r25*t46*t55*zh*(-1.0_r64/2.0_r64)
            pcs(22) = ipi*rhop*t22*t34*r25*zh*(t39+rho*t4*6.0_r64-rho*t6*3.0_r64-rho*t17*6.0_r64-t3*t4*4.0_r64+t3*t6*4.0_r64 &
                    +chi*t2*t17*2.0_r64+chi*t6*t19)*(-1.0_r64/8.0_r64)
            pcs(23) = ipi*nrp*rho*r15*zh*cmplx(0.0_r64,0.75_r64,r64)
            pcs(24) = ipi*t33*r15*t56*zh*(-1.0_r64/4.0_r64)
            pcs(25) = ipi*rhop*t15*r25*t46*t55*(-1.0_r64/8.0_r64)
            pcs(26) = ipi*(rhop*t15*r25*t46*t55)/2.0_r64
            pcs(27) = ipi*(rhop*t15*t22*t34*r25*(t16+chi*t6*4.0_r64+chi*t17-nrp*t3*4.0_r64))/8.0_r64
          end block
          do md = 0, M
            rn = real(md,r64); n2 = rn*rn; ifn2 = 1.0_r64/(4.0_r64*n2-1.0_r64)
            vk = vka(md); ve = vea(md)
            coK(1)=pcs(1)+n2*pcs(2);   coK(3)=pcs(8)+n2*pcs(9);   coK(5)=pcs(14)+n2*pcs(15)
            coK(7)=pcs(20)+n2*pcs(21); coK(9)=pcs(25)+n2*pcs(26)
            coK(2)=rn*pcs(5); coK(4)=rn*pcs(11); coK(6)=rn*pcs(18); coK(8)=rn*pcs(23)
            coE(1)=ifn2*(pcs(3)+n2*pcs(4)); coE(5)=ifn2*(pcs(16)+n2*pcs(17))
            coE(2)=rn*ifn2*(pcs(6)+n2*pcs(7)); coE(4)=rn*ifn2*(pcs(12)+n2*pcs(13))
            coE(3)=pcs(10); coE(7)=pcs(22); coE(9)=pcs(27)
            coE(6)=rn*pcs(19); coE(8)=rn*pcs(24)
            ! 3x3 block entries  Be_** = (coK*vk + coE*ve)*ws   (row r/t/z = 0/1/2*nt, col r/t/z = 0/1/2*nso)
            Be_rr(jp) = (coK(1)*vk+coE(1)*ve)*ws;  Be_rt(jp) = (coK(2)*vk+coE(2)*ve)*ws;  Be_rz(jp) = (coK(3)*vk+coE(3)*ve)*ws
            Be_tr(jp) = (coK(4)*vk+coE(4)*ve)*ws;  Be_tt(jp) = (coK(5)*vk+coE(5)*ve)*ws;  Be_tz(jp) = (coK(6)*vk+coE(6)*ve)*ws
            Be_zr(jp) = (coK(7)*vk+coE(7)*ve)*ws;  Be_zt(jp) = (coK(8)*vk+coE(8)*ve)*ws;  Be_zz(jp) = (coK(9)*vk+coE(9)*ve)*ws
            A(0*nt+i, 0*nso+jj, md+1) = Be_rr(jp);  A(0*nt+i, 1*nso+jj, md+1) = Be_rt(jp);  A(0*nt+i, 2*nso+jj, md+1) = Be_rz(jp)
            A(1*nt+i, 0*nso+jj, md+1) = Be_tr(jp);  A(1*nt+i, 1*nso+jj, md+1) = Be_tt(jp);  A(1*nt+i, 2*nso+jj, md+1) = Be_tz(jp)
            A(2*nt+i, 0*nso+jj, md+1) = Be_zr(jp);  A(2*nt+i, 1*nso+jj, md+1) = Be_zt(jp);  A(2*nt+i, 2*nso+jj, md+1) = Be_zz(jp)
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
              select case (e)                                  ! fold Gpc into the (row,col) 3x3 block of A
              case (1); do ia=1,nkc; do l=1,p; A(0*nt+nidx(ci(ia)), 0*nso+cols+l, md+1) = A(0*nt+nidx(ci(ia)), 0*nso+cols+l, md+1) + Gpc(ia,l); end do; end do  ! Be_rr
              case (2); do ia=1,nkc; do l=1,p; A(0*nt+nidx(ci(ia)), 1*nso+cols+l, md+1) = A(0*nt+nidx(ci(ia)), 1*nso+cols+l, md+1) + Gpc(ia,l); end do; end do  ! Be_rt
              case (3); do ia=1,nkc; do l=1,p; A(0*nt+nidx(ci(ia)), 2*nso+cols+l, md+1) = A(0*nt+nidx(ci(ia)), 2*nso+cols+l, md+1) + Gpc(ia,l); end do; end do  ! Be_rz
              case (4); do ia=1,nkc; do l=1,p; A(1*nt+nidx(ci(ia)), 0*nso+cols+l, md+1) = A(1*nt+nidx(ci(ia)), 0*nso+cols+l, md+1) + Gpc(ia,l); end do; end do  ! Be_tr
              case (5); do ia=1,nkc; do l=1,p; A(1*nt+nidx(ci(ia)), 1*nso+cols+l, md+1) = A(1*nt+nidx(ci(ia)), 1*nso+cols+l, md+1) + Gpc(ia,l); end do; end do  ! Be_tt
              case (6); do ia=1,nkc; do l=1,p; A(1*nt+nidx(ci(ia)), 2*nso+cols+l, md+1) = A(1*nt+nidx(ci(ia)), 2*nso+cols+l, md+1) + Gpc(ia,l); end do; end do  ! Be_tz
              case (7); do ia=1,nkc; do l=1,p; A(2*nt+nidx(ci(ia)), 0*nso+cols+l, md+1) = A(2*nt+nidx(ci(ia)), 0*nso+cols+l, md+1) + Gpc(ia,l); end do; end do  ! Be_zr
              case (8); do ia=1,nkc; do l=1,p; A(2*nt+nidx(ci(ia)), 1*nso+cols+l, md+1) = A(2*nt+nidx(ci(ia)), 1*nso+cols+l, md+1) + Gpc(ia,l); end do; end do  ! Be_zt
              case (9); do ia=1,nkc; do l=1,p; A(2*nt+nidx(ci(ia)), 2*nso+cols+l, md+1) = A(2*nt+nidx(ci(ia)), 2*nso+cols+l, md+1) + Gpc(ia,l); end do; end do  ! Be_zz
              end select
            end do
          end do
        end if
        ! ---- inner far: naive dlp9far on the FINER subpanel nodes, fold via Lc ----
        do i = 1, nk
          if (ikc(i)) cycle
          rt = real(znear(i),r64); zti = aimag(znear(i))
          do jp = 1, p                                          ! carrier_all + md-indep pieces ONCE per node
            rhop = real(Ypb(jp),r64); zh = zti - aimag(Ypb(jp)); rho = rt
            nrp = real(-ic*dYp(jp)/abs(dYp(jp)),r64); nzp = aimag(-ic*dYp(jp)/abs(dYp(jp)))
            rr2 = (rho-rhop)**2 + zh*zh; chi = 1.0_r64 + rr2/(2.0_r64*rho*rhop)
            call modal_green_all_far_r64(chi, M, vka, vea)
            vkmat(jp,0:M) = vka(0:M); vemat(jp,0:M) = vea(0:M)
            block
              real(r64) :: ipi, r15, r25
              real(r64) :: t2,t3,t4,t5,t6,t7,t8,t11,t12,t13,t15,t16,t17,t18,t19,t20,t21,t22,t23,t26,t27,t28,t29,t30,t31,t32
              real(r64) :: t33,t34,t36,t39,t40,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55
              complex(r64) :: t24,t25,t35,t37,t41,t56
              ipi = 1.0_r64/acos(-1.0_r64)
              t2=chi*rho; t3=chi*rhop; t4=nrp*rhop; t5=rho*rhop; t6=nzp*zh; t7=chi+1.0_r64; t8=chi*chi
              t11=rho*rho; t12=t11*rho; t13=rhop*rhop; t15=zh*zh; t16=nrp*rho*3.0_r64; t18=chi-1.0_r64
              t17=nrp*t2; t19=-t2; t20=-t3; t21=-t4; t22=1.0_r64/t7; t23=t8-1.0_r64
              r15=1.0_r64/t5**1.5_r64; r25=1.0_r64/t5**2.5_r64
              t24=t4*ic; t25=t6*ic; t26=nrp*t11*3.0_r64; t27=t4*t13*3.0_r64; t31=t6*t11*3.0_r64; t32=t6*t13*3.0_r64
              t33=1.0_r64/t18; t40=rhop*t2*t6*10.0_r64; t41=nrp*t11*3.0_r64*ic; t45=rhop*t2*t4*16.0_r64
              t48=t2*t2*t6; t49=t3*t3*t6; t50=chi*t2*t3*t6*2.0_r64
              t28=t11*t17*6.0_r64; t29=rho+t20; t30=rhop+t19; t34=t33*t33; t35=-t24; t36=-t27; t37=t17*ic
              t39=t3*t17*2.0_r64; t44=t2*t2*t17*2.0_r64; t46=1.0_r64/t23; t47=-t40; t52=t3*t4*t20; t53=t3*t3*t17*4.0_r64
              t55=t6+t17+t21; t51=-t44; t54=-t53; t56=t25+t35+t37
              pcsm(1,jp) = ipi*rhop*r25*t46*(nrp*t12*3.0_r64+t4*t5*7.0_r64-t5*t6*4.0_r64-rho*t2*t4*11.0_r64+rho*t2*t6-rho*t2*t17*2.0_r64 &
                      +rhop*t3*t6+rhop*t4*t20-t2*t3*t4*4.0_r64+t2*t3*t6*2.0_r64+t2*t3*t17*8.0_r64)*(-1.0_r64/8.0_r64)
              pcsm(2,jp) = ipi*rhop*r25*t46*t55*(t5+rho*t19+rhop*t20+t2*t3)*(-1.0_r64/2.0_r64)
              pcsm(3,jp) = ipi*rhop*t22*t34*r25*(t28+t31+t32+t36+t45+t47+t48+t49+t50+t51+t52+t54-t4*t11*3.0_r64 &
                      -t2*t2*t4*17.0_r64+chi*t2*t3*t17*8.0_r64)*(-1.0_r64/8.0_r64)
              pcsm(4,jp) = ipi*(rhop*t22*t34*r25*(t28+t31+t32+t36+t45+t47+t48+t49+t50+t51+t52+t54-t4*t11*6.0_r64 &
                      -t2*t2*t4*11.0_r64+chi*t2*t3*t17*5.0_r64))/2.0_r64
              pcsm(5,jp) = ipi*(r15*(t41+rhop*t4*3.0_r64*ic-rhop*t6*3.0_r64*ic-t2*t4*6.0_r64*ic))/4.0_r64
              pcsm(6,jp) = ipi*(t33*r15*(rho*t4*(-4.0_r64*ic)+rho*t25+rho*t37-t3*t6*ic+t3*t17*2.0_r64*ic+t3*t24))/4.0_r64
              pcsm(7,jp) = ipi*(-t29*t33*r15*t56)
              pcsm(8,jp) = ipi*rhop*r25*t46*zh*(t39-rho*t4*4.0_r64+rho*t6+rho*t17+t3*t4+t6*t20)*(-1.0_r64/8.0_r64)
              pcsm(9,jp) = ipi*(rhop*t29*r25*t46*t55*zh)/2.0_r64
              pcsm(10,jp) = ipi*(rhop*t22*t34*r25*zh*(t26+chi*t39+rhop*t4*3.0_r64-rhop*t6*3.0_r64-t2*t4*10.0_r64+t2*t6*4.0_r64 &
                      +t2*t17+nrp*t3*t3+chi*t6*t20))/8.0_r64
              pcsm(11,jp) = ipi*rhop*r15*(t4*(-2.0_r64*ic)+t17*2.0_r64*ic+t25)*(3.0_r64/4.0_r64)
              pcsm(12,jp) = ipi*(t13*t33*r25*(t41-rhop*t6*ic+rhop*t24-t2*t4*2.0_r64*ic-t2*t17*2.0_r64*ic+t2*t25))/4.0_r64
              pcsm(13,jp) = ipi*t13*t30*t33*r25*t56
              pcsm(14,jp) = ipi*rhop*r15*(t17*3.0_r64+t55)*(-1.0_r64/4.0_r64)
              pcsm(15,jp) = ipi*rhop*r15*t55*(-1.0_r64/2.0_r64)
              pcsm(16,jp) = ipi*(rhop*t33*r15*(t16-chi*t6-chi*t17*4.0_r64+nrp*t3))/4.0_r64
              pcsm(17,jp) = ipi*rhop*t33*r15*(t16-chi*t6*2.0_r64-chi*t17*5.0_r64+nrp*t3*2.0_r64)*(-1.0_r64/2.0_r64)
              pcsm(18,jp) = ipi*t4*r15*zh*cmplx(0.0_r64,0.75_r64,r64)
              pcsm(19,jp) = ipi*t13*t33*r25*t56*zh*(-1.0_r64/4.0_r64)
              pcsm(20,jp) = ipi*rhop*r25*t46*zh*(t26+rhop*t4-rhop*t6-t2*t4*2.0_r64+t2*t6-t2*t17*2.0_r64)*(-1.0_r64/8.0_r64)
              pcsm(21,jp) = ipi*rhop*t30*r25*t46*t55*zh*(-1.0_r64/2.0_r64)
              pcsm(22,jp) = ipi*rhop*t22*t34*r25*zh*(t39+rho*t4*6.0_r64-rho*t6*3.0_r64-rho*t17*6.0_r64-t3*t4*4.0_r64+t3*t6*4.0_r64 &
                      +chi*t2*t17*2.0_r64+chi*t6*t19)*(-1.0_r64/8.0_r64)
              pcsm(23,jp) = ipi*nrp*rho*r15*zh*cmplx(0.0_r64,0.75_r64,r64)
              pcsm(24,jp) = ipi*t33*r15*t56*zh*(-1.0_r64/4.0_r64)
              pcsm(25,jp) = ipi*rhop*t15*r25*t46*t55*(-1.0_r64/8.0_r64)
              pcsm(26,jp) = ipi*(rhop*t15*r25*t46*t55)/2.0_r64
              pcsm(27,jp) = ipi*(rhop*t15*t22*t34*r25*(t16+chi*t6*4.0_r64+chi*t17-nrp*t3*4.0_r64))/8.0_r64
            end block
          end do
          do md = 0, M
            rn = real(md,r64); n2 = rn*rn; ifn2 = 1.0_r64/(4.0_r64*n2-1.0_r64)
            do jp = 1, p
              vk = vkmat(jp,md); ve = vemat(jp,md)
              coK(1)=pcsm(1,jp)+n2*pcsm(2,jp);   coK(3)=pcsm(8,jp)+n2*pcsm(9,jp);   coK(5)=pcsm(14,jp)+n2*pcsm(15,jp)
              coK(7)=pcsm(20,jp)+n2*pcsm(21,jp); coK(9)=pcsm(25,jp)+n2*pcsm(26,jp)
              coK(2)=rn*pcsm(5,jp); coK(4)=rn*pcsm(11,jp); coK(6)=rn*pcsm(18,jp); coK(8)=rn*pcsm(23,jp)
              coE(1)=ifn2*(pcsm(3,jp)+n2*pcsm(4,jp)); coE(5)=ifn2*(pcsm(16,jp)+n2*pcsm(17,jp))
              coE(2)=rn*ifn2*(pcsm(6,jp)+n2*pcsm(7,jp)); coE(4)=rn*ifn2*(pcsm(12,jp)+n2*pcsm(13,jp))
              coE(3)=pcsm(10,jp); coE(7)=pcsm(22,jp); coE(9)=pcsm(27,jp)
              coE(6)=rn*pcsm(19,jp); coE(8)=rn*pcsm(24,jp)
              ! 3x3 block entries  Be_** = (coK*vk + coE*ve)*wsp   (per subpanel node jp, folded below)
              Be_rr(jp) = (coK(1)*vk+coE(1)*ve)*wsp(jp);  Be_rt(jp) = (coK(2)*vk+coE(2)*ve)*wsp(jp);  Be_rz(jp) = (coK(3)*vk+coE(3)*ve)*wsp(jp)
              Be_tr(jp) = (coK(4)*vk+coE(4)*ve)*wsp(jp);  Be_tt(jp) = (coK(5)*vk+coE(5)*ve)*wsp(jp);  Be_tz(jp) = (coK(6)*vk+coE(6)*ve)*wsp(jp)
              Be_zr(jp) = (coK(7)*vk+coE(7)*ve)*wsp(jp);  Be_zt(jp) = (coK(8)*vk+coE(8)*ve)*wsp(jp);  Be_zz(jp) = (coK(9)*vk+coE(9)*ve)*wsp(jp)
            end do
            do l = 1, p                                        ! fold each named block (sum over subpanel nodes via Lc) into A
              A(0*nt+nidx(i), 0*nso+cols+l, md+1) = A(0*nt+nidx(i), 0*nso+cols+l, md+1) + sum(Be_rr(1:p)*Lc(:,l))
              A(0*nt+nidx(i), 1*nso+cols+l, md+1) = A(0*nt+nidx(i), 1*nso+cols+l, md+1) + sum(Be_rt(1:p)*Lc(:,l))
              A(0*nt+nidx(i), 2*nso+cols+l, md+1) = A(0*nt+nidx(i), 2*nso+cols+l, md+1) + sum(Be_rz(1:p)*Lc(:,l))
              A(1*nt+nidx(i), 0*nso+cols+l, md+1) = A(1*nt+nidx(i), 0*nso+cols+l, md+1) + sum(Be_tr(1:p)*Lc(:,l))
              A(1*nt+nidx(i), 1*nso+cols+l, md+1) = A(1*nt+nidx(i), 1*nso+cols+l, md+1) + sum(Be_tt(1:p)*Lc(:,l))
              A(1*nt+nidx(i), 2*nso+cols+l, md+1) = A(1*nt+nidx(i), 2*nso+cols+l, md+1) + sum(Be_tz(1:p)*Lc(:,l))
              A(2*nt+nidx(i), 0*nso+cols+l, md+1) = A(2*nt+nidx(i), 0*nso+cols+l, md+1) + sum(Be_zr(1:p)*Lc(:,l))
              A(2*nt+nidx(i), 1*nso+cols+l, md+1) = A(2*nt+nidx(i), 1*nso+cols+l, md+1) + sum(Be_zt(1:p)*Lc(:,l))
              A(2*nt+nidx(i), 2*nso+cols+l, md+1) = A(2*nt+nidx(i), 2*nso+cols+l, md+1) + sum(Be_zz(1:p)*Lc(:,l))
            end do
          end do
        end do
      end do
    end do
    deallocate(tin, joa, ci, zc, nidx, znear, ikq, ikc, C1, C2, C3, C4, As, Ad, A1, A2, A3, A4, Gcq, Gpc)
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
                call modal_green_r64(chi, 0_8, vk, ve, Fn, An, dFn)
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
    real(r64)  :: vka(0:M), vea(0:M), Fna(0:M), Ana(0:M), dFna(0:M), vkmat(p,0:M), vemat(p,0:M)
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
    complex(r64) :: coK(9), coE(9), pcs(36), pcsm(36,p)
    complex(r64) :: Be_rr(p), Be_rt(p), Be_rz(p), &       ! 9 entries of the 3x3 block, per panel node
                    Be_tr(p), Be_tt(p), Be_tz(p), &
                    Be_zr(p), Be_zt(p), Be_zz(p)
    real(r64)  :: n2, ifn2, n4, nr, nz
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
          call modal_green_all_far_r64(chi, M, vka, vea)        ! vk/ve only (D' far recipe)
          nr = nrt; nz = nzt                                    ! target normal feeds the D' coK/coE pieces
          block
            real(r64) :: ipi, t2, t3, t4, t5, t6, t7, t8, t9, t10, t15, t16, t18, t19, t22, t27, t28, t20, t25, t26, t29, t30, t31, t32
            real(r64) :: t33, t34, t36, t37, t41, t42, t43, t44, t45, t46, t47, t48, t49, t50, t51, t52, t54, t55, t56, t57, t58, t60
            real(r64) :: t69, t72, t73, t74, t75, t76, t77, t78, t79, t80, t82, t83, t84, t87, t92, t93, t109, t110, t111, t112, t114
            real(r64) :: t115, t116, t117, t118, t119, t125, t128, t129, t134, t136, t137, t138, t139, t141, t144, t145, t147, t185
            real(r64) :: t188, t190, t191, t192, t198, t201, t205, t208, t210, t211, t212, t213, t214, t215, t231, t250, t252, t253
            real(r64) :: t254, t257, t259, t260, t262, t263, t265, t266, t274, t278, t279, t280, t281, t291, t292, t35, t38, t39, t40
            real(r64) :: t53, t59, t61, t65, t66, t67, t70, t71, t81, t85, t91, t102, t108, t113, t120, t121, t122, t123, t124, t127
            real(r64) :: t130, t131, t132, t133, t135, t140, t142, t146, t157, t159, t162, t176, t177, t178, t179, t180, t181, t182
            real(r64) :: t183, t184, t187, t189, t193, t203, t204, t206, t207, t229, t230, t236, t241, t245, t247, t248, t249, t251
            real(r64) :: t258, t261, t264, t272, t273, t275, t277, t288, t289, t290, t294, t295, t296, t297, t298, t299, t301, t68
            real(r64) :: t126, t143, t158, t165, t170, t171, t186, t202, t244, t255, t256, t276, t293, t300, et1, et2, et3, et4, et5
            real(r64) :: et6, et7, et8, et9, et10
            complex(r64) :: t62, t63, t64, t94, t95, t96, t97, t98, t99, t100, t101, t103, t104, t105, t106, t107, t160, t161, t167, t169
            complex(r64) :: t172, t173, t175, t220, t221, t232, t233, t234, t235, t148, t149, t150, t151, t152, t153, t154, t155, t156
            complex(r64) :: t163, t164, t166, t168, t174, t194, t195, t216, t217, t218, t219, t222, t223, t224, t226, t227, t228, t237
            complex(r64) :: t238, t239, t240, t242, t269, t270, t282, t283, t284, t285, t196, t197, t225, t243, t267, t268, t271, t286
            complex(r64) :: t287
            ipi = 1.0_r64/acos(-1.0_r64)
            t2 = chi*rho
            t3 = chi*rhop
            t4 = nr*rho
            t5 = nrp*rhop
            t6 = rho*rhop
            t7 = nz*zh
            t8 = nzp*zh
            t9 = chi+1.0_r64
            t10 = chi**2
            t15 = rho**2
            t16 = rho**3
            t18 = rhop**2
            t19 = rhop**3
            t22 = zh**2
            t27 = chi-1.0_r64
            t28 = 1.0_r64/rho
            t20 = t18**2
            t25 = nr*t3
            t26 = nrp*t2
            t29 = 1.0_r64/t15
            t30 = -t2
            t31 = -t3
            t32 = -t5
            t33 = 1.0_r64/t9
            t34 = t10-1.0_r64
            t36 = nr*t2*t8
            t37 = nrp*t3*t7
            t41 = 1.0_r64/t27**2
            t42 = 1.0_r64/t27**3
            t43 = chi*t7*t8
            t44 = nz*nzp*rhop*t2*2.0_r64
            t45 = nrp*t4*t16*1.5e+1_r64
            t46 = nr*t5*t19*1.5e+1_r64
            t47 = nz*nzp*t6*t18*3.0_r64
            t48 = nz*nzp*t6*t15*3.0_r64
            t49 = nr*t5*t19*zh*5.0_r64
            t50 = nrp*rho*t7*t22*5.0_r64
            t51 = nr*rhop*t8*t22*5.0_r64
            t52 = t5*t7*t19*5.0_r64
            t54 = 1.0_r64/t6**(3.0_r64/2.0_r64)
            t55 = 1.0_r64/t6**(5.0_r64/2.0_r64)
            t56 = 1.0_r64/t6**(7.0_r64/2.0_r64)
            t57 = nr*t2*t5*t18
            t58 = rho*t2*t4*t5
            t60 = nz*t2*t5*t19
            t62 = nrp*t4*t15*(5.0_r64*ic)
            t63 = nr*t5*t18*(5.0_r64*ic)
            t64 = t7*t8*zh*(1.0_r64*ic)
            t69 = nrp*t7*t16*1.5e+1_r64
            t72 = nr*t8*t19*1.5e+1_r64
            t73 = chi*nz*nzp*t2*t3*2.0_r64
            t74 = t4*t5*t15*zh*5.0_r64
            t75 = t4*t6*t8*2.0_r64
            t76 = t5*t6*t7*2.0_r64
            t77 = rho*t4*t6*t8*5.0_r64
            t78 = nzp*t6*t7*t18*2.0_r64
            t79 = t6*t7*t8*4.0_r64
            t80 = t6*t7*t8*zh*5.0_r64
            t82 = nr*t2**3*t5
            t83 = nzp*rhop*t2*t4*t6
            t84 = nz*rho*t2*t5*t6
            t87 = nz*nzp*t2*t6**2
            t92 = rho*t2*t5*t7*2.0_r64
            t93 = nzp*t2*t7*t19*2.0_r64
            t94 = nz*nzp*rhop*t6*(1.0_r64*ic)
            t95 = nz*nzp*rho*t6*(1.0_r64*ic)
            t96 = t4*t8*zh*(1.0_r64*ic)
            t97 = rho*t4*t8*(4.0_r64*ic)
            t98 = nrp*t7*t15*(5.0_r64*ic)
            t99 = t5*t7*zh*(1.0_r64*ic)
            t100 = rhop*t5*t7*(4.0_r64*ic)
            t101 = nr*t8*t18*(5.0_r64*ic)
            t103 = rho*t7*t8*(4.0_r64*ic)
            t104 = rhop*t7*t8*(4.0_r64*ic)
            t105 = t4*t5*zh*(1.0e+1_r64*ic)
            t106 = rhop*t4*t8*(8.0_r64*ic)
            t107 = rho*t5*t7*(8.0_r64*ic)
            t109 = nzp*t4*t6*t18*4.0_r64
            t110 = nz*t5*t6**2*4.0_r64
            t111 = nz*nzp*rhop*t6**2*4.0_r64
            t112 = t2*t4*t5*t15*1.0e+1_r64
            t114 = nrp*rho*t4*t22*5.0_r64
            t115 = nr*rhop*t5*t22*5.0_r64
            t116 = nr*t8*t19*zh*5.0_r64
            t117 = t5*t7*t18*zh*5.0_r64
            t118 = t7*t8*t18*zh*5.0_r64
            t119 = t7*t8*t19*5.0_r64
            t125 = t2**3*t4*t5*2.0_r64
            t128 = nz*nzp*t2*t3*t18*2.0_r64
            t129 = nz*nzp*t2**2*t6*2.0_r64
            t134 = rhop*t3**2*t5*t7*3.0_r64
            t136 = t4*t6*t8*zh*5.0_r64
            t137 = t5*t7*t15*zh*5.0_r64
            t138 = nzp*t6**2*t7*2.0_r64
            t139 = rho*t6*t7*t8*5.0_r64
            t141 = nr*t2*t5*t22*1.7e+1_r64
            t144 = t2*t5*t7*t15*1.0e+1_r64
            t145 = t2*t5*t7*t18*2.6e+1_r64
            t147 = rhop*t2*t7*t8*zh*8.0_r64
            t160 = nz*nzp*t2*t18*(1.0_r64*ic)
            t161 = nz*nzp*t2*t6*(1.0_r64*ic)
            t167 = t2*t4*t8*(4.0_r64*ic)
            t169 = t3*t5*t7*(4.0_r64*ic)
            t172 = t2*t7*t8*(4.0_r64*ic)
            t173 = t3*t7*t8*(4.0_r64*ic)
            t175 = t2*t5*t7*(1.2e+1_r64*ic)
            t185 = nz*nzp*rhop*t2*t6*1.0e+1_r64
            t188 = rhop*t3*t5*t7*zh*8.0_r64
            t190 = rhop*t3*t7*t8*zh*8.0_r64
            t191 = t4*t8*t18*zh*1.0e+1_r64
            t192 = t5*t6*t7*zh*1.0e+1_r64
            t198 = nr*nzp*t2**3*t18
            t201 = t2**2*t4*t5*zh*3.0_r64
            t205 = t2*t3*t5*t7*2.0_r64
            t208 = rhop*t2**2*t4*t8*3.0_r64
            t210 = t2**3*t5*t7*2.0_r64
            t211 = nzp*rhop*t2*t6*t7*2.0_r64
            t212 = nzp*t2*t3*t7*t18*2.0_r64
            t213 = t2*t3*t7*t8*4.0_r64
            t214 = t2*t3*t7*t8*zh*3.0_r64
            t215 = nzp*rhop*t2*t3**2*t7*2.0_r64
            t220 = nr*rhop*t2*t5*(2.0e+1_r64*ic)
            t221 = t2*t4*t5*(2.0e+1_r64*ic)
            t231 = t2*t3**2*t7*t8
            t232 = nz*nzp*rhop*t2*t3*(1.0_r64*ic)
            t233 = nz*nzp*rhop*t2**2*(1.0_r64*ic)
            t234 = nz*nzp*t2*t3**2*(1.0_r64*ic)
            t235 = nz*nzp*t2**2*t3*(1.0_r64*ic)
            t250 = t3**2*t5*t7*zh*3.0_r64
            t252 = t3**2*t7*t8*zh*3.0_r64
            t253 = rhop*t3**2*t7*t8*3.0_r64
            t254 = rhop*t2*t4*t5*zh*2.6e+1_r64
            t257 = t2*t3**2*t5*t7*(-1.0_r64)*2.0_r64
            t259 = rhop*t2*t4*t8*zh*8.0_r64
            t260 = rho*t2*t5*t7*zh*1.0e+1_r64
            t262 = rhop*t2*t5*t7*zh*1.7e+1_r64
            t263 = t2*t4*t8*t18*2.6e+1_r64
            t265 = t2*t7*t8*t18*1.7e+1_r64
            t266 = t2*t3*t7*t8*4.4e+1_r64
            t274 = t2**2*t5*t7*zh*3.0_r64
            t278 = t2*t3*t5*t7*zh*6.0_r64
            t279 = nzp*t2**2*t7*t18*2.0_r64
            t280 = rhop*t2**2*t7*t8*3.0_r64
            t281 = nzp*rhop*t2**2*t3*t7*2.0_r64
            t291 = chi*nz*nzp*t2**2*t3**2*(-1.0_r64)*2.0_r64
            t292 = nz*nzp*rhop*t2**2*t3*1.2e+1_r64
            t35 = t33**2
            t38 = -t25
            t39 = rho+t31
            t40 = rhop+t30
            t53 = nrp*t7*t31
            t59 = nr*nzp*t2*t20
            t61 = nz*nzp*t2*t20
            t65 = t5*t19*t25*1.0e+1_r64
            t66 = t22*t43*8.0_r64
            t67 = 1.0_r64/t34
            t70 = -t51
            t71 = -t52
            t81 = t2*t3*t5*t25
            t85 = nzp*t2*t3*t18*t25
            t91 = t18*t36*2.0_r64
            t102 = -t72
            t108 = rhop*t3**2*t5*t25*2.0_r64
            t113 = -t73
            t120 = t5*t18*t25*zh*1.0e+1_r64
            t121 = t22*t36*8.0_r64
            t122 = t22*t37*8.0_r64
            t123 = t8*t19*t25*1.0e+1_r64
            t124 = -t76
            t127 = t82*4.0_r64
            t130 = rhop*t3*t5*t25*zh*3.0_r64
            t131 = t3**2*t5*t25*zh*2.0_r64
            t132 = chi*t7*t22*t26*3.0_r64
            t133 = chi*t8*t22*t25*3.0_r64
            t135 = rhop*t3**2*t8*t25*2.0_r64
            t140 = t58*zh*1.0e+1_r64
            t142 = t57*zh*2.6e+1_r64
            t146 = -t93
            t148 = rho*t4*t26*(8.0_r64*ic)
            t149 = rhop*t5*t25*(8.0_r64*ic)
            t150 = -t97
            t151 = -t99
            t152 = -t101
            t153 = -t103
            t154 = -t104
            t155 = -t105
            t156 = -t106
            t157 = chi*t22*t25*t26
            t159 = -t111
            t162 = -t118
            t163 = t4*t26*zh*(1.0_r64*ic)
            t164 = t5*t25*zh*(1.0_r64*ic)
            t166 = t7*t26*zh*(1.0_r64*ic)
            t168 = t8*t25*zh*(1.0_r64*ic)
            t174 = rhop*t36*(1.2e+1_r64*ic)
            t176 = t3*t5*t25*t30
            t177 = nr*t2**2*t5*t30
            t178 = -t125
            t179 = nz*rhop*t3**2*t5*t30
            t180 = -t128
            t181 = -t129
            t182 = nz*nzp*t3**2*t18*t30
            t183 = nz*nzp*t3**3*t30
            t184 = nz*nzp*t2**2*t3*t30
            t187 = -t134
            t189 = t8*t18*t25*zh*1.0e+1_r64
            t193 = -t141
            t194 = chi*t2*t8*t25*(4.0_r64*ic)
            t195 = chi*t3*t7*t26*(4.0_r64*ic)
            t203 = t82*zh*2.0_r64
            t204 = rhop*t2*t36*2.0_r64
            t206 = t2*t3*t8*t25*2.0_r64
            t207 = t2*t3*t7*t26*2.0_r64
            t216 = t2*t4*t26*(1.0_r64*ic)
            t217 = nr*t2**2*t26*(4.0_r64*ic)
            t218 = t3*t5*t25*(1.0_r64*ic)
            t219 = nrp*t3**2*t25*(4.0_r64*ic)
            t222 = -t160
            t223 = -t161
            t224 = rho*t7*t26*(8.0_r64*ic)
            t226 = -t169
            t227 = rhop*t8*t25*(8.0_r64*ic)
            t228 = -t175
            t229 = t2*t3*t8*t25*zh
            t230 = t3**2*t7*t26*zh
            t236 = -t185
            t237 = t2*t7*t26*(1.0_r64*ic)
            t238 = chi*t2*t7*t26*(4.0_r64*ic)
            t239 = t3*t8*t25*(1.0_r64*ic)
            t240 = chi*t3*t8*t25*(4.0_r64*ic)
            t241 = -t191
            t242 = t25*t26*zh*(8.0_r64*ic)
            t245 = t8+t26+t32
            t247 = nr*t2*t22*t26*3.0_r64
            t248 = nrp*t3*t22*t25*3.0_r64
            t249 = rhop*t3*t8*t25*zh*3.0_r64
            t251 = t3**2*t8*t25*zh*2.0_r64
            t258 = -t210
            t261 = t18*t36*zh*1.7e+1_r64
            t264 = -t213
            t269 = -t220
            t270 = -t221
            t272 = t2**2*t5*t25*zh*2.0_r64
            t273 = rhop*t2*t36*zh*3.0_r64
            t275 = rhop*t2**2*t8*t25*2.0_r64
            t277 = rhop*t2*t8*t25*zh*6.0_r64
            t282 = t3*t25*t26*(8.0_r64*ic)
            t283 = t2*t25*t26*(8.0_r64*ic)
            t284 = -t232
            t285 = -t233
            t288 = rhop*t2**2*t5*t25*1.2e+1_r64
            t289 = nz*rhop*t2**2*t5*t30
            t290 = nz*nzp*t2**2*t18*t30
            t294 = -t252
            t295 = -t254
            t296 = -t262
            t297 = -t263
            t298 = -t265
            t299 = -t266
            t301 = -t279
            t68 = t67**2
            t126 = t81*4.0_r64
            t143 = -t91
            t158 = -t108
            t165 = -t120
            t170 = -t122
            t171 = -t123
            t186 = -t133
            t196 = -t148
            t197 = -t149
            t202 = t81*zh*2.0_r64
            t225 = -t168
            t243 = -t194
            t244 = t4+t7+t38
            t255 = -t204
            t256 = -t207
            t267 = -t216
            t268 = -t218
            t271 = -t224
            t276 = t207*zh
            t286 = -t237
            t287 = -t240
            t293 = -t251
            t300 = -t277
            et1 = t45+t46+t47+t48-t57*1.24e+2_r64-t58*1.24e+2_r64+t69+t102+t126+t127+t180+t181+t183+t184+t236+t291+t292+t299+t18*t36*9.1e+1_r64+t2**2*t36*3.0_r64-t3**2*t37*3.0_r64+t2**2*t43*3.0_r64+t3**2*t43*3.0_r64-nr*t2**3*t26*2.0_r64-nrp*t3**3*t25*2.0_r64-rhop*t2*t36*6.5e+1_r64+t4*t5*t6*5.0e+1_r64-t4*t6*t8*3.3e+1_r64+t5*t6*t7*3.3e+1_r64-t6*t7*t8*2.1e+1_r64+rho*t2*t4*t8*2.9e+1_r64-rho*t2*t5*t7*9.1e+1_r64+rho*t2*t7*t8*2.9e+1_r64+rho*t2*t4*t26*1.9e+1_r64+rho*t2*t7*t26*1.9e+1_r64-rhop*t3*t5*t7*2.9e+1_r64+rhop*t3*t7*t8*2.9e+1_r64+rhop*t3*t5*t25*1.9e+1_r64-rhop*t3*t8*t25*1.9e+1_r64+t2*t3*t5*t7*6.5e+1_r64+t2*t3*t8*t25-t2*t3*t25*t26*5.2e+1_r64
            et2 = t3*t7*t26*t30-chi*t2**2*t7*t26*2.0_r64+chi*t2**2*t8*t25*2.0_r64-chi*t3**2*t7*t26*2.0_r64+chi*t3**2*t8*t25*2.0_r64-chi*t2**2*t25*t26*8.0_r64-chi*t3**2*t25*t26*8.0_r64+nr*rhop*t2**2*t5*1.78e+2_r64+t2*t3*t7*t8*t10-t2*t3*t7*t10*t26*4.0_r64+t2*t3*t8*t10*t25*4.0_r64+t2*t3*t10*t25*t26*1.6e+1_r64
            et3 = t45+t46+t47+t48-t57*1.26e+2_r64-t58*1.26e+2_r64+t69+t102+t126+t127+t180+t181+t183+t184+t236+t291+t292+t299+t18*t36*9.3e+1_r64+t2**2*t36*2.0_r64-t3**2*t37*2.0_r64+t2**2*t43*2.0_r64+t3**2*t43*2.0_r64-nr*t2**3*t26*3.0_r64-nrp*t3**3*t25*3.0_r64-rhop*t2*t36*6.6e+1_r64+t4*t5*t6*5.4e+1_r64-t4*t6*t8*3.4e+1_r64+t5*t6*t7*3.4e+1_r64-t6*t7*t8*2.2e+1_r64+rho*t2*t4*t8*3.0e+1_r64-rho*t2*t5*t7*9.3e+1_r64+rho*t2*t7*t8*3.0e+1_r64+rho*t2*t4*t26*2.0e+1_r64+rho*t2*t7*t26*2.0e+1_r64-rhop*t3*t5*t7*3.0e+1_r64+rhop*t3*t7*t8*3.0e+1_r64+rhop*t3*t5*t25*2.0e+1_r64-rhop*t3*t8*t25*2.0e+1_r64+t2*t3*t5*t7*6.6e+1_r64-t2*t3*t25*t26*4.6e+1_r64
            et4 = chi*t2**2*t7*t26*(-1.0_r64)*3.0_r64+chi*t2**2*t8*t25*4.0_r64-chi*t3**2*t7*t26*4.0_r64+chi*t3**2*t8*t25*3.0_r64-chi*t2**2*t25*t26*6.0_r64-chi*t3**2*t25*t26*6.0_r64+nr*rhop*t2**2*t5*1.72e+2_r64+t2*t3*t7*t8*t10*2.0_r64-t2*t3*t7*t10*t26*3.0_r64+t2*t3*t8*t10*t25*3.0_r64+t2*t3*t10*t25*t26*1.2e+1_r64
            et5 = -t109-t229*6.0_r64+chi*t230+t57*zh*9.1e+1_r64+t58*zh*2.9e+1_r64+t72*zh+t81*zh+t82*zh*3.0_r64+t3**2*t37*zh*3.0_r64-t3**2*t43*zh*3.0_r64-t18*t36*zh*5.8e+1_r64+nz*t5*t6*t15*3.0_r64+nz*t5*t6*t18*3.0_r64-nzp*t6*t7*t18*6.0_r64-nr*t5*t19*zh*1.5e+1_r64+rhop*t2*t36*zh*2.3e+1_r64-t4*t5*t6*zh*3.3e+1_r64+t4*t6*t8*zh*9.0_r64-t5*t6*t7*zh*2.1e+1_r64+t6*t7*t8*zh*9.0_r64+nr*nzp*t2**2*t19*8.0_r64-nz*t2**2*t5*t6*2.0_r64-nz*t2*t3**3*t26*2.0_r64+nz*t3**3*t5*t30+nzp*t2*t3**3*t7*2.0_r64+nrp*t3**3*t25*zh*2.0_r64-nr*rhop*t2**2*t5*zh*6.5e+1_r64-nz*rhop*t2*t5*t6*1.0e+1_r64+nzp*rhop*t2*t6*t7*8.0_r64
            et6 = nz*t2*t3*t5*t18*(-1.0_r64)*2.0_r64+nzp*t2*t3*t7*t18*4.0_r64+rho*t2*t5*t7*zh*2.9e+1_r64+rhop*t3*t5*t7*zh*2.9e+1_r64-rhop*t3*t7*t8*zh*2.9e+1_r64-rhop*t3*t5*t25*zh*1.9e+1_r64+rhop*t3*t8*t25*zh*1.9e+1_r64-t2*t3*t5*t7*zh*4.4e+1_r64+t2*t3*t7*t8*zh*2.3e+1_r64+t2*t3*t7*t26*zh*3.0_r64+t2*t3*t25*t26*zh*2.0_r64-chi*t3**2*t8*t25*zh*2.0_r64+chi*t3**2*t25*t26*zh*4.0_r64+nz*rhop*t2**2*t3*t5*1.2e+1_r64-nzp*rhop*t2**2*t3*t7*8.0_r64-nzp*rhop*t2**2*t3*t25*4.0_r64+nz*t2**2*t3*t5*t30
            et7 = t110+t230*6.0_r64+chi*t276-nzp*t4*t6**2*3.0_r64-nzp*t6**2*t7*6.0_r64-nzp*t4*t20*3.0_r64+nzp*t2*t4*t19*1.0e+1_r64+nzp*t2*t7*t19*8.0_r64+nzp*t2*t19*t25*2.0_r64-t4*t5*t15*zh*1.5e+1_r64-t4*t5*t18*zh*3.3e+1_r64-t5*t7*t15*zh*1.5e+1_r64+t4*t8*t18*zh*2.1e+1_r64-t5*t7*t18*zh*9.0_r64+t7*t8*t18*zh*9.0_r64-t2*t3*t43*zh*3.0_r64+t5*t18*t25*zh*2.9e+1_r64-t8*t18*t25*zh*2.9e+1_r64+nzp*rhop*t2**3*t25-nz*t2**2*t5*t18*8.0_r64+nzp*t2**2*t4*t18*2.0_r64+nzp*t2**2*t7*t18*4.0_r64-nzp*t2**2*t18*t25*1.2e+1_r64-t2**2*t4*t5*zh*1.9e+1_r64-t2**2*t5*t7*zh*1.9e+1_r64-t3**2*t5*t7*zh*2.3e+1_r64+t3**2*t7*t8*zh*2.3e+1_r64
            et8 = t2**2*t5*t25*zh+t3**2*t5*t25*zh*3.0_r64-t2**2*t8*t25*zh*3.0_r64-t3**2*t8*t25*zh*3.0_r64+t2**2*t25*t26*zh*2.0_r64+t3**2*t25*t26*zh*2.0_r64+nz*t2**2*t3**2*t5*4.0_r64+nzp*t2**2*t3**2*t7*2.0_r64+nzp*t2**2*t3**2*t25*2.0_r64+rhop*t2*t4*t5*zh*9.1e+1_r64-rhop*t2*t4*t8*zh*2.9e+1_r64+rhop*t2*t5*t7*zh*5.8e+1_r64-rhop*t2*t7*t8*zh*2.9e+1_r64-rhop*t2*t5*t25*zh*6.5e+1_r64+rhop*t2*t8*t25*zh*4.4e+1_r64-nzp*rhop*t2*t3**2*t7*8.0_r64+nzp*rhop*t2*t3**2*t25+chi*t2*t3*t25*t26*zh*4.0_r64+chi*t3*t8*t25*t30*zh
            et9 = t159-t231*1.6e+1_r64+t288+t3*t132+t3*t157+t3*t176+t3*t206+t3*t248+t4*t5*t19*3.0_r64-t4*t8*t19*6.0_r64-t3*t22*t37*2.3e+1_r64+t3*t22*t43*2.3e+1_r64+rho*t4*t5*t6*3.0_r64+rho*t5*t6*t7*6.0_r64-rhop*t4*t5*t22*2.1e+1_r64+rhop*t4*t8*t22*9.0_r64-rhop*t5*t7*t22*9.0_r64+rhop*t7*t8*t22*9.0_r64+rhop*t5*t22*t25*2.9e+1_r64-rhop*t8*t22*t25*2.9e+1_r64-t2*t4*t5*t18*1.0e+1_r64+t2*t4*t8*t18*8.0_r64-t2*t5*t7*t18*8.0_r64+t2*t4*t5*t22*2.9e+1_r64+t2*t7*t8*t18*1.6e+1_r64+t2*t5*t7*t22*2.9e+1_r64-t2*t5*t18*t25*2.0_r64+t2*t8*t18*t25*4.0_r64-t2*t5*t22*t25*4.4e+1_r64+t2*t8*t22*t25*2.3e+1_r64+t2*t22*t25*t26*3.0_r64+nz*nzp*t2**2*t19*8.0_r64-rhop*t2**2*t4*t5*2.0_r64
            et10 = rhop*t2**2*t5*t7*(-1.0_r64)*4.0_r64-rhop*t2**2*t8*t25*8.0_r64+t2*t3**2*t5*t7*8.0_r64-t2*t3**2*t7*t26*2.0_r64-t2*t3**2*t25*t26*2.0_r64+t2**2*t5*t25*t30-nz*nzp*rhop*t2**2*t3**2*4.0_r64-chi*t3*t8*t22*t25*3.0_r64
            pcs(1) = ipi*((t56*t68*(t61+t65+t71+t77+t87+t112+t119+t135+t139+t144+t145+t158+t159+t171+t178+t182+t187+t208+t231+t253+t257+t258+t275+t280+t290+t297+t298-t3*t81*8.0_r64-t2**3*t5*t25*8.0_r64-t4*t5*t19*2.8e+1_r64+t4*t8*t19*1.9e+1_r64-rho*t4*t5*t6*2.8e+1_r64-rho*t5*t6*t7*1.9e+1_r64+t2*t4*t5*t18*6.8e+1_r64+t2*t5*t18*t25*4.0_r64+t2*t8*t18*t25+nz*nzp*t2**2*t19*6.0_r64+rhop*t2**2*t4*t5*4.0_r64-rhop*t2**2*t5*t25*3.6e+1_r64-t2*t3**2*t7*t26*4.0_r64+t2*t3**2*t8*t25*4.0_r64+t2*t3**2*t25*t26*1.6e+1_r64-nz*nzp*rhop*t2**2*t3**2*2.0_r64+rhop*t2*t5*t7*t30))/8.0_r64)
            pcs(2) = ipi*(t56*t68*(t61+t65+t71+t77+t87+t112+t119+t135+t139+t144+t145+t158+t171+t178+t182+t187+t208+t231+t253+t257+t258+t275+t280+t288+t290+t297+t298+t3*t81*1.0e+1_r64+t2**3*t5*t25*1.0e+1_r64-t4*t5*t19*1.0e+1_r64+t4*t8*t19*1.0e+1_r64-rho*t4*t5*t6*1.0e+1_r64-rho*t5*t6*t7*1.0e+1_r64+t2*t4*t5*t18*4.4e+1_r64-t2*t5*t18*t25*3.2e+1_r64+t2*t8*t18*t25*1.9e+1_r64-nz*nzp*rhop*t6**2-rhop*t2**2*t4*t5*3.2e+1_r64-rhop*t2**2*t5*t7*1.9e+1_r64+t2*t3**2*t7*t26*5.0_r64-t2*t3**2*t8*t25*5.0_r64-t2*t3**2*t25*t26*8.0_r64+nz*nzp*rhop*t2**2*t3**2)*(-1.0_r64/2.0_r64))
            pcs(3) = ipi*((t28*t35*t42*t55*(et1+et2))/8.0_r64)
            pcs(4) = ipi*(t28*t35*t42*t55*(et3+et4)*(-1.0_r64/2.0_r64))
            pcs(5) = ipi*(t28*t33*t41*t55*t244*t245*(t6+rho*t30+rhop*t31+t2*t3)*2.0_r64)
            pcs(6) = ipi*(t55*t67*(rhop*t36*(-1.0_r64)*(2.0_r64*ic)+nr*t2**2*t5*(1.6e+1_r64*ic)-nr*t5*t18*(9.0_r64*ic)+nr*t8*t18*(9.0_r64*ic)-rho*t4*t5*(1.9e+1_r64*ic)+rho*t4*t8*(1.0_r64*ic)-rho*t5*t7*(1.0e+1_r64*ic)+rho*t7*t8*(1.0_r64*ic)+rho*t4*t26*(1.0_r64*ic)+rho*t7*t26*(1.0_r64*ic)+t3*t5*t7*(1.0_r64*ic)-t3*t7*t8*(1.0_r64*ic)+t3*t5*t25*(8.0_r64*ic)+t3*t7*t26*(8.0_r64*ic)-t3*t8*t25*(8.0_r64*ic)-t3*t25*t26*(2.0e+1_r64*ic)-nz*nzp*rhop*t6*(3.0_r64*ic)+nr*rhop*t2*t5*(2.3e+1_r64*ic)+nz*nzp*rhop*t2*t3*(3.0_r64*ic))*(-1.0_r64/4.0_r64))
            pcs(7) = ipi*(t39*t55*t67*t244*t245*(1.0_r64*ic))
            pcs(8) = ipi*(t33*t41*t55*(t62+t95+t98+t100+t154+t156+t167+t172+t195+t197+t219+t222+t227+t228+t234+t267+t270+t283+t285+t286+t287+rhop*t4*t5*(1.1e+1_r64*ic)+t2*t5*t25*(5.0_r64*ic)-chi*t3*t25*t26*(4.0_r64*ic))*(-1.0_r64/4.0_r64))
            pcs(9) = ipi*(t33*t41*t55*(t62+t95+t98+t100+t154+t156+t167+t172+t195+t197+t219+t222+t227+t228+t234+t267+t270+t283+t285+t286+t287+rhop*t4*t5*(8.0_r64*ic)+t2*t5*t25*(1.1e+1_r64*ic)-chi*t3*t25*t26*(7.0_r64*ic)))
            pcs(10) = ipi*((t56*t68*(t60+t74+t84-t110+t117+t131+t137+t138+t146+t147+t162+t165+t179+t189+t201+t215+t230+t241+t250+t259+t272+t274+t289+t293+t294+t295+t296+t300+t301+t4*t5*t18*zh*1.9e+1_r64+nz*t2**2*t5*t18*6.0_r64+t3**2*t25*t26*zh*4.0_r64-nz*t2**2*t3**2*t5*2.0_r64+rhop*t2*t5*t25*zh))/8.0_r64)
            pcs(11) = ipi*(t56*t68*(t60+t74+t84+t117+t131+t137+t138+t146+t147+t162+t165+t179+t189+t201+t215+t230+t241+t250+t259+t272+t274+t289+t293+t294+t295+t296+t300+t301+nz*t6**2*t32+t4*t5*t18*zh*1.0e+1_r64-t3**2*t25*t26*zh*5.0_r64+nz*t2**2*t3**2*t5+rhop*t2*t5*t25*zh*1.9e+1_r64)*(-1.0_r64/2.0_r64))
            pcs(12) = ipi*(t35*t42*t56*(et5+et6)*(-1.0_r64/8.0_r64))
            pcs(13) = ipi*(rhop*t33*t39*t41*t56*t244*t245*zh*(-1.0_r64/2.0_r64))
            pcs(14) = ipi*(t18*t56*t67*(nrp*t4*t15*(9.0_r64*ic)+nrp*t7*t15*(9.0_r64*ic)+rhop*t4*t5*(1.9e+1_r64*ic)-rhop*t4*t8*(1.0e+1_r64*ic)+rhop*t5*t7*(1.0_r64*ic)-rhop*t7*t8*(1.0_r64*ic)-rhop*t5*t25*(1.0_r64*ic)+rhop*t8*t25*(1.0_r64*ic)-t2*t4*t5*(2.3e+1_r64*ic)+t2*t4*t8*(1.0_r64*ic)-t2*t5*t7*(2.0_r64*ic)+t2*t7*t8*(1.0_r64*ic)-t2*t4*t26*(8.0_r64*ic)-t2*t5*t25*(1.6e+1_r64*ic)-t2*t7*t26*(8.0_r64*ic)+t2*t8*t25*(8.0_r64*ic)+t2*t25*t26*(2.0e+1_r64*ic)+nz*nzp*rho*t6*(3.0_r64*ic)-nz*nzp*rhop*t2**2*(3.0_r64*ic))*(-1.0_r64/4.0_r64))
            pcs(15) = ipi*(t18*t40*t56*t67*t244*t245*(-1.0_r64)*(1.0_r64*ic))
            pcs(16) = ipi*((t29*t33*t41*t54*(t63+t94+t107+t150+t152+t153+t173+t174+t196+t217+t223+t226+t235+t238+t239+t243+t268+t269+t271+t282+t284+nr*t2**2*t5*(5.0_r64*ic)+rho*t4*t5*(1.1e+1_r64*ic)-chi*t2*t25*t26*(4.0_r64*ic)))/4.0_r64)
            pcs(17) = ipi*(-t29*t33*t41*t54*(t63+t94+t107+t150+t152+t153+t173+t174+t196+t217+t223+t226+t235+t238+t239+t243+t268+t269+t271+t282+t284+nr*t2**2*t5*(1.1e+1_r64*ic)+rho*t4*t5*(8.0_r64*ic)-chi*t2*t25*t26*(7.0_r64*ic)))
            pcs(18) = ipi*((t28*t54*t67*(t36+t43+t53-chi*t7*t26*4.0_r64+chi*t8*t25*4.0_r64+chi*t25*t26*1.6e+1_r64+nz*nzp*t6*2.0_r64+nr*rhop*t5*5.0_r64-nr*rhop*t8*5.0_r64+nrp*rho*t4*5.0_r64+nrp*rho*t7*5.0_r64-nr*t2*t5*1.8e+1_r64-nr*t2*t26*4.0_r64-nrp*t3*t25*4.0_r64-nz*nzp*t2*t3*2.0_r64))/8.0_r64)
            pcs(19) = ipi*(t28*t54*t67*(t36+t43+t53+chi*t7*t26*5.0_r64-chi*t8*t25*5.0_r64-chi*t25*t26*8.0_r64-nz*nzp*t6-nr*rhop*t5*4.0_r64+nr*rhop*t8*4.0_r64-nrp*rho*t4*4.0_r64-nrp*rho*t7*4.0_r64+nr*t2*t5*6.0_r64+nr*t2*t26*5.0_r64+nrp*t3*t25*5.0_r64+nz*nzp*t2*t3)*(-1.0_r64/2.0_r64))
            pcs(20) = ipi*((rhop*t33*t41*t55*(t44+t113+chi*t36+chi*t53+t4*t5*6.0_r64+t4*t8*3.0_r64-t5*t7*3.0_r64+t7*t8*3.0_r64+t4*t26*8.0_r64+t5*t25*8.0_r64+t7*t26*8.0_r64-t8*t25*8.0_r64-t25*t26*3.0e+1_r64+t7*t8*t10-t7*t10*t26*4.0_r64+t8*t10*t25*4.0_r64+t10*t25*t26*1.6e+1_r64-chi*nr*t2*t26*4.0_r64-chi*nrp*t3*t25*4.0_r64))/8.0_r64)
            pcs(21) = ipi*(rhop*t33*t41*t55*(t44+t113+chi*t36*2.0_r64-chi*t37*2.0_r64+t4*t5*4.0_r64+t4*t8*2.0_r64-t5*t7*2.0_r64+t7*t8*2.0_r64+t4*t26*7.0_r64+t5*t25*7.0_r64+t7*t26*7.0_r64-t8*t25*7.0_r64-t25*t26*2.4e+1_r64+t7*t8*t10*2.0_r64-t7*t10*t26*3.0_r64+t8*t10*t25*3.0_r64+t10*t25*t26*1.2e+1_r64-chi*nr*t2*t26*3.0_r64-chi*nrp*t3*t25*3.0_r64)*(-1.0_r64/2.0_r64))
            pcs(22) = ipi*((rhop*t55*t244*t245*2.0_r64)/t27)
            pcs(23) = ipi*(t29*t54*t67*(t64+t96+t151+t155+t163+t164+t166+t225+t242-nz*t2**2*t5*(3.0_r64*ic)+nz*t5*t15*(3.0_r64*ic))*(-1.0_r64/4.0_r64))
            pcs(24) = ipi*(t29*t54*t67*t244*t245*zh*(1.0_r64*ic))
            pcs(25) = ipi*((t28*t33*t41*t55*(t3*t168+t63*zh+t173*zh+rhop*t36*zh*(4.0_r64*ic)-nz*rhop*t5*t6*(1.0_r64*ic)+nzp*rhop*t6*t7*(2.0_r64*ic)+nz*t2*t5*t6*(1.0_r64*ic)-nr*t8*t18*zh*(5.0_r64*ic)+rho*t4*t5*zh*(4.0_r64*ic)+rho*t5*t7*zh*(4.0_r64*ic)-t3*t5*t7*zh*(4.0_r64*ic)-t3*t5*t25*zh*(1.0_r64*ic)+t3*t25*t26*zh*(4.0_r64*ic)-nz*t2**2*t3*t5*(1.0_r64*ic)+nz*rhop*t2*t3*t5*(1.0_r64*ic)-nzp*rhop*t2*t3*t7*(2.0_r64*ic)-nr*rhop*t2*t5*zh*(1.2e+1_r64*ic)))/4.0_r64)
            pcs(26) = ipi*(t56*t68*(t49+t78-t80+t85+t109-t116+t130-t136-t140-t142-t188+t190+t192+t198+t202+t203-t211-t212-t214-t249-t260+t261-t273+t276+t278+t281+nr*nzp*t20*t30+t4*t5*t6*zh*1.9e+1_r64-nr*nzp*t2**2*t19*6.0_r64+nr*rhop*t2**2*t5*zh+nzp*rhop*t4*t6*t30+t2*t3*t25*t26*zh*4.0_r64+t3*t8*t25*t30*zh+nzp*rhop*t2**2*t3*t25*2.0_r64)*(-1.0_r64/8.0_r64))
            pcs(27) = ipi*(t56*t68*(-t49+t59-t78+t80+t83+t116-t130+t136+t140+t142+t188-t190-t192-t202-t203+t211+t212+t214+t229+t249+t260-t261+t273-t278-t281-nzp*t4*t6*t18-t4*t5*t6*zh*1.0e+1_r64-nr*rhop*t2**2*t5*zh*1.9e+1_r64+nzp*t3*t18*t25*t30-t2*t3*t7*t26*zh*2.0_r64+t2*t3*t25*t26*zh*5.0_r64+nr*nzp*t2**2*t18*t30+nzp*rhop*t2**2*t3*t25)*(-1.0_r64/2.0_r64))
            pcs(28) = ipi*((t35*t42*t56*(et7+et8))/8.0_r64)
            pcs(29) = ipi*((rhop*t33*t40*t41*t56*t244*t245*zh)/2.0_r64)
            pcs(30) = ipi*(t55*t67*(t64+t96+t151+t155+t163+t164+t166+t225+t242-nzp*t4*t18*(3.0_r64*ic)+nzp*rhop*t2*t25*(3.0_r64*ic))*(-1.0_r64/4.0_r64))
            pcs(31) = ipi*(t55*t67*t244*t245*zh*(1.0_r64*ic))
            pcs(32) = ipi*((t33*t41*t55*(t36*zh*(4.0_r64*ic)-t37*zh*(4.0_r64*ic)+t43*zh*(4.0_r64*ic)+nzp*t4*t6*(1.0_r64*ic)+nzp*t6*t7*(2.0_r64*ic)-nr*nzp*t2*t18*(1.0_r64*ic)-chi*t7*t26*zh*(1.0_r64*ic)+chi*t25*t26*zh*(4.0_r64*ic)-nzp*t2*t3*t7*(2.0_r64*ic)+nzp*t2*t3*t25*(1.0_r64*ic)+nr*rhop*t5*zh*(4.0_r64*ic)-nr*rhop*t8*zh*(4.0_r64*ic)+nrp*rho*t4*zh*(5.0_r64*ic)+nrp*rho*t7*zh*(5.0_r64*ic)-nr*t2*t5*zh*(1.2e+1_r64*ic)-nr*t2*t26*zh*(1.0_r64*ic)-nr*nzp*rhop*t2**2*(1.0_r64*ic)))/4.0_r64)
            pcs(33) = ipi*((t28*t55*t68*(t50+t57+t58+t66+t70+t75+t79+t92+t114+t115+t121+t124+t132+t143+t157+t170+t176+t177+t186+t193+t205+t206+t247+t248+t255+t256+t264-t4*t5*t6*4.0_r64-t2*t3*t25*t26*2.0_r64+nr*rhop*t2**2*t5*6.0_r64))/8.0_r64)
            pcs(34) = ipi*(t28*t55*t68*(t50+t57+t58+t66+t70+t75+t79+t92+t114+t115+t121+t124+t132+t143+t157+t170+t176+t177+t186+t193+t205+t206+t247+t248+t255+t256+t264+t4*t6*t32+t2*t3*t25*t26)*(-1.0_r64/2.0_r64))
            pcs(35) = ipi*(t35*t42*t56*(et9+et10)*(-1.0_r64/8.0_r64))
            pcs(36) = ipi*(rhop*t22*t33*t41*t56*t244*t245*(-1.0_r64/2.0_r64))
          end block
          do md = 0, M
            rn = real(md,r64); n2 = rn*rn; ifn2 = 1.0_r64/(4.0_r64*n2-1.0_r64); n4 = n2*n2
            vk = vka(md); ve = vea(md)
            coK(1)=pcs(1)+n2*pcs(2);   coK(3)=pcs(10)+n2*pcs(11); coK(5)=pcs(18)+n2*pcs(19)
            coK(7)=pcs(26)+n2*pcs(27); coK(9)=pcs(33)+n2*pcs(34)
            coK(2)=rn*(pcs(6)+n2*pcs(7));   coK(4)=rn*(pcs(14)+n2*pcs(15))
            coK(6)=rn*(pcs(23)+n2*pcs(24)); coK(8)=rn*(pcs(30)+n2*pcs(31))
            coE(1)=ifn2*(pcs(3)+n2*pcs(4)+n4*pcs(5)); coE(5)=ifn2*(pcs(20)+n2*pcs(21)+n4*pcs(22))
            coE(3)=pcs(12)+n2*pcs(13); coE(7)=pcs(28)+n2*pcs(29); coE(9)=pcs(35)+n2*pcs(36)
            coE(2)=rn*ifn2*(pcs(8)+n2*pcs(9)); coE(4)=rn*ifn2*(pcs(16)+n2*pcs(17))
            coE(6)=rn*pcs(25); coE(8)=rn*pcs(32)
            ! 3x3 block entries  Be_** = (coK*vk + coE*ve)*ws
            Be_rr(jp) = (coK(1)*vk+coE(1)*ve)*ws;  Be_rt(jp) = (coK(2)*vk+coE(2)*ve)*ws;  Be_rz(jp) = (coK(3)*vk+coE(3)*ve)*ws
            Be_tr(jp) = (coK(4)*vk+coE(4)*ve)*ws;  Be_tt(jp) = (coK(5)*vk+coE(5)*ve)*ws;  Be_tz(jp) = (coK(6)*vk+coE(6)*ve)*ws
            Be_zr(jp) = (coK(7)*vk+coE(7)*ve)*ws;  Be_zt(jp) = (coK(8)*vk+coE(8)*ve)*ws;  Be_zz(jp) = (coK(9)*vk+coE(9)*ve)*ws
            A(0*nt+i, 0*nso+jj, md+1) = Be_rr(jp);  A(0*nt+i, 1*nso+jj, md+1) = Be_rt(jp);  A(0*nt+i, 2*nso+jj, md+1) = Be_rz(jp)
            A(1*nt+i, 0*nso+jj, md+1) = Be_tr(jp);  A(1*nt+i, 1*nso+jj, md+1) = Be_tt(jp);  A(1*nt+i, 2*nso+jj, md+1) = Be_tz(jp)
            A(2*nt+i, 0*nso+jj, md+1) = Be_zr(jp);  A(2*nt+i, 1*nso+jj, md+1) = Be_zt(jp);  A(2*nt+i, 2*nso+jj, md+1) = Be_zz(jp)
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
          do jp = 1, p                                          ! carrier_all + md-indep D' pieces ONCE per node
            rhop = real(Ypb(jp),r64); zh = zti - aimag(Ypb(jp)); rho = rt
            nrp = real(-ic*dYp(jp)/abs(dYp(jp)),r64); nzp = aimag(-ic*dYp(jp)/abs(dYp(jp)))
            nr = nrt; nz = nzt
            rr2 = (rho-rhop)**2 + zh*zh; chi = 1.0_r64 + rr2/(2.0_r64*rho*rhop)
            call modal_green_all_far_r64(chi, M, vka, vea)
            vkmat(jp,0:M) = vka(0:M); vemat(jp,0:M) = vea(0:M)
            block
              real(r64) :: ipi, t2, t3, t4, t5, t6, t7, t8, t9, t10, t15, t16, t18, t19, t22, t27, t28, t20, t25, t26, t29, t30, t31, t32
              real(r64) :: t33, t34, t36, t37, t41, t42, t43, t44, t45, t46, t47, t48, t49, t50, t51, t52, t54, t55, t56, t57, t58, t60
              real(r64) :: t69, t72, t73, t74, t75, t76, t77, t78, t79, t80, t82, t83, t84, t87, t92, t93, t109, t110, t111, t112, t114
              real(r64) :: t115, t116, t117, t118, t119, t125, t128, t129, t134, t136, t137, t138, t139, t141, t144, t145, t147, t185
              real(r64) :: t188, t190, t191, t192, t198, t201, t205, t208, t210, t211, t212, t213, t214, t215, t231, t250, t252, t253
              real(r64) :: t254, t257, t259, t260, t262, t263, t265, t266, t274, t278, t279, t280, t281, t291, t292, t35, t38, t39, t40
              real(r64) :: t53, t59, t61, t65, t66, t67, t70, t71, t81, t85, t91, t102, t108, t113, t120, t121, t122, t123, t124, t127
              real(r64) :: t130, t131, t132, t133, t135, t140, t142, t146, t157, t159, t162, t176, t177, t178, t179, t180, t181, t182
              real(r64) :: t183, t184, t187, t189, t193, t203, t204, t206, t207, t229, t230, t236, t241, t245, t247, t248, t249, t251
              real(r64) :: t258, t261, t264, t272, t273, t275, t277, t288, t289, t290, t294, t295, t296, t297, t298, t299, t301, t68
              real(r64) :: t126, t143, t158, t165, t170, t171, t186, t202, t244, t255, t256, t276, t293, t300, et1, et2, et3, et4, et5
              real(r64) :: et6, et7, et8, et9, et10
              complex(r64) :: t62, t63, t64, t94, t95, t96, t97, t98, t99, t100, t101, t103, t104, t105, t106, t107, t160, t161, t167, t169
              complex(r64) :: t172, t173, t175, t220, t221, t232, t233, t234, t235, t148, t149, t150, t151, t152, t153, t154, t155, t156
              complex(r64) :: t163, t164, t166, t168, t174, t194, t195, t216, t217, t218, t219, t222, t223, t224, t226, t227, t228, t237
              complex(r64) :: t238, t239, t240, t242, t269, t270, t282, t283, t284, t285, t196, t197, t225, t243, t267, t268, t271, t286
              complex(r64) :: t287
              ipi = 1.0_r64/acos(-1.0_r64)
              t2 = chi*rho
              t3 = chi*rhop
              t4 = nr*rho
              t5 = nrp*rhop
              t6 = rho*rhop
              t7 = nz*zh
              t8 = nzp*zh
              t9 = chi+1.0_r64
              t10 = chi**2
              t15 = rho**2
              t16 = rho**3
              t18 = rhop**2
              t19 = rhop**3
              t22 = zh**2
              t27 = chi-1.0_r64
              t28 = 1.0_r64/rho
              t20 = t18**2
              t25 = nr*t3
              t26 = nrp*t2
              t29 = 1.0_r64/t15
              t30 = -t2
              t31 = -t3
              t32 = -t5
              t33 = 1.0_r64/t9
              t34 = t10-1.0_r64
              t36 = nr*t2*t8
              t37 = nrp*t3*t7
              t41 = 1.0_r64/t27**2
              t42 = 1.0_r64/t27**3
              t43 = chi*t7*t8
              t44 = nz*nzp*rhop*t2*2.0_r64
              t45 = nrp*t4*t16*1.5e+1_r64
              t46 = nr*t5*t19*1.5e+1_r64
              t47 = nz*nzp*t6*t18*3.0_r64
              t48 = nz*nzp*t6*t15*3.0_r64
              t49 = nr*t5*t19*zh*5.0_r64
              t50 = nrp*rho*t7*t22*5.0_r64
              t51 = nr*rhop*t8*t22*5.0_r64
              t52 = t5*t7*t19*5.0_r64
              t54 = 1.0_r64/t6**(3.0_r64/2.0_r64)
              t55 = 1.0_r64/t6**(5.0_r64/2.0_r64)
              t56 = 1.0_r64/t6**(7.0_r64/2.0_r64)
              t57 = nr*t2*t5*t18
              t58 = rho*t2*t4*t5
              t60 = nz*t2*t5*t19
              t62 = nrp*t4*t15*(5.0_r64*ic)
              t63 = nr*t5*t18*(5.0_r64*ic)
              t64 = t7*t8*zh*(1.0_r64*ic)
              t69 = nrp*t7*t16*1.5e+1_r64
              t72 = nr*t8*t19*1.5e+1_r64
              t73 = chi*nz*nzp*t2*t3*2.0_r64
              t74 = t4*t5*t15*zh*5.0_r64
              t75 = t4*t6*t8*2.0_r64
              t76 = t5*t6*t7*2.0_r64
              t77 = rho*t4*t6*t8*5.0_r64
              t78 = nzp*t6*t7*t18*2.0_r64
              t79 = t6*t7*t8*4.0_r64
              t80 = t6*t7*t8*zh*5.0_r64
              t82 = nr*t2**3*t5
              t83 = nzp*rhop*t2*t4*t6
              t84 = nz*rho*t2*t5*t6
              t87 = nz*nzp*t2*t6**2
              t92 = rho*t2*t5*t7*2.0_r64
              t93 = nzp*t2*t7*t19*2.0_r64
              t94 = nz*nzp*rhop*t6*(1.0_r64*ic)
              t95 = nz*nzp*rho*t6*(1.0_r64*ic)
              t96 = t4*t8*zh*(1.0_r64*ic)
              t97 = rho*t4*t8*(4.0_r64*ic)
              t98 = nrp*t7*t15*(5.0_r64*ic)
              t99 = t5*t7*zh*(1.0_r64*ic)
              t100 = rhop*t5*t7*(4.0_r64*ic)
              t101 = nr*t8*t18*(5.0_r64*ic)
              t103 = rho*t7*t8*(4.0_r64*ic)
              t104 = rhop*t7*t8*(4.0_r64*ic)
              t105 = t4*t5*zh*(1.0e+1_r64*ic)
              t106 = rhop*t4*t8*(8.0_r64*ic)
              t107 = rho*t5*t7*(8.0_r64*ic)
              t109 = nzp*t4*t6*t18*4.0_r64
              t110 = nz*t5*t6**2*4.0_r64
              t111 = nz*nzp*rhop*t6**2*4.0_r64
              t112 = t2*t4*t5*t15*1.0e+1_r64
              t114 = nrp*rho*t4*t22*5.0_r64
              t115 = nr*rhop*t5*t22*5.0_r64
              t116 = nr*t8*t19*zh*5.0_r64
              t117 = t5*t7*t18*zh*5.0_r64
              t118 = t7*t8*t18*zh*5.0_r64
              t119 = t7*t8*t19*5.0_r64
              t125 = t2**3*t4*t5*2.0_r64
              t128 = nz*nzp*t2*t3*t18*2.0_r64
              t129 = nz*nzp*t2**2*t6*2.0_r64
              t134 = rhop*t3**2*t5*t7*3.0_r64
              t136 = t4*t6*t8*zh*5.0_r64
              t137 = t5*t7*t15*zh*5.0_r64
              t138 = nzp*t6**2*t7*2.0_r64
              t139 = rho*t6*t7*t8*5.0_r64
              t141 = nr*t2*t5*t22*1.7e+1_r64
              t144 = t2*t5*t7*t15*1.0e+1_r64
              t145 = t2*t5*t7*t18*2.6e+1_r64
              t147 = rhop*t2*t7*t8*zh*8.0_r64
              t160 = nz*nzp*t2*t18*(1.0_r64*ic)
              t161 = nz*nzp*t2*t6*(1.0_r64*ic)
              t167 = t2*t4*t8*(4.0_r64*ic)
              t169 = t3*t5*t7*(4.0_r64*ic)
              t172 = t2*t7*t8*(4.0_r64*ic)
              t173 = t3*t7*t8*(4.0_r64*ic)
              t175 = t2*t5*t7*(1.2e+1_r64*ic)
              t185 = nz*nzp*rhop*t2*t6*1.0e+1_r64
              t188 = rhop*t3*t5*t7*zh*8.0_r64
              t190 = rhop*t3*t7*t8*zh*8.0_r64
              t191 = t4*t8*t18*zh*1.0e+1_r64
              t192 = t5*t6*t7*zh*1.0e+1_r64
              t198 = nr*nzp*t2**3*t18
              t201 = t2**2*t4*t5*zh*3.0_r64
              t205 = t2*t3*t5*t7*2.0_r64
              t208 = rhop*t2**2*t4*t8*3.0_r64
              t210 = t2**3*t5*t7*2.0_r64
              t211 = nzp*rhop*t2*t6*t7*2.0_r64
              t212 = nzp*t2*t3*t7*t18*2.0_r64
              t213 = t2*t3*t7*t8*4.0_r64
              t214 = t2*t3*t7*t8*zh*3.0_r64
              t215 = nzp*rhop*t2*t3**2*t7*2.0_r64
              t220 = nr*rhop*t2*t5*(2.0e+1_r64*ic)
              t221 = t2*t4*t5*(2.0e+1_r64*ic)
              t231 = t2*t3**2*t7*t8
              t232 = nz*nzp*rhop*t2*t3*(1.0_r64*ic)
              t233 = nz*nzp*rhop*t2**2*(1.0_r64*ic)
              t234 = nz*nzp*t2*t3**2*(1.0_r64*ic)
              t235 = nz*nzp*t2**2*t3*(1.0_r64*ic)
              t250 = t3**2*t5*t7*zh*3.0_r64
              t252 = t3**2*t7*t8*zh*3.0_r64
              t253 = rhop*t3**2*t7*t8*3.0_r64
              t254 = rhop*t2*t4*t5*zh*2.6e+1_r64
              t257 = t2*t3**2*t5*t7*(-1.0_r64)*2.0_r64
              t259 = rhop*t2*t4*t8*zh*8.0_r64
              t260 = rho*t2*t5*t7*zh*1.0e+1_r64
              t262 = rhop*t2*t5*t7*zh*1.7e+1_r64
              t263 = t2*t4*t8*t18*2.6e+1_r64
              t265 = t2*t7*t8*t18*1.7e+1_r64
              t266 = t2*t3*t7*t8*4.4e+1_r64
              t274 = t2**2*t5*t7*zh*3.0_r64
              t278 = t2*t3*t5*t7*zh*6.0_r64
              t279 = nzp*t2**2*t7*t18*2.0_r64
              t280 = rhop*t2**2*t7*t8*3.0_r64
              t281 = nzp*rhop*t2**2*t3*t7*2.0_r64
              t291 = chi*nz*nzp*t2**2*t3**2*(-1.0_r64)*2.0_r64
              t292 = nz*nzp*rhop*t2**2*t3*1.2e+1_r64
              t35 = t33**2
              t38 = -t25
              t39 = rho+t31
              t40 = rhop+t30
              t53 = nrp*t7*t31
              t59 = nr*nzp*t2*t20
              t61 = nz*nzp*t2*t20
              t65 = t5*t19*t25*1.0e+1_r64
              t66 = t22*t43*8.0_r64
              t67 = 1.0_r64/t34
              t70 = -t51
              t71 = -t52
              t81 = t2*t3*t5*t25
              t85 = nzp*t2*t3*t18*t25
              t91 = t18*t36*2.0_r64
              t102 = -t72
              t108 = rhop*t3**2*t5*t25*2.0_r64
              t113 = -t73
              t120 = t5*t18*t25*zh*1.0e+1_r64
              t121 = t22*t36*8.0_r64
              t122 = t22*t37*8.0_r64
              t123 = t8*t19*t25*1.0e+1_r64
              t124 = -t76
              t127 = t82*4.0_r64
              t130 = rhop*t3*t5*t25*zh*3.0_r64
              t131 = t3**2*t5*t25*zh*2.0_r64
              t132 = chi*t7*t22*t26*3.0_r64
              t133 = chi*t8*t22*t25*3.0_r64
              t135 = rhop*t3**2*t8*t25*2.0_r64
              t140 = t58*zh*1.0e+1_r64
              t142 = t57*zh*2.6e+1_r64
              t146 = -t93
              t148 = rho*t4*t26*(8.0_r64*ic)
              t149 = rhop*t5*t25*(8.0_r64*ic)
              t150 = -t97
              t151 = -t99
              t152 = -t101
              t153 = -t103
              t154 = -t104
              t155 = -t105
              t156 = -t106
              t157 = chi*t22*t25*t26
              t159 = -t111
              t162 = -t118
              t163 = t4*t26*zh*(1.0_r64*ic)
              t164 = t5*t25*zh*(1.0_r64*ic)
              t166 = t7*t26*zh*(1.0_r64*ic)
              t168 = t8*t25*zh*(1.0_r64*ic)
              t174 = rhop*t36*(1.2e+1_r64*ic)
              t176 = t3*t5*t25*t30
              t177 = nr*t2**2*t5*t30
              t178 = -t125
              t179 = nz*rhop*t3**2*t5*t30
              t180 = -t128
              t181 = -t129
              t182 = nz*nzp*t3**2*t18*t30
              t183 = nz*nzp*t3**3*t30
              t184 = nz*nzp*t2**2*t3*t30
              t187 = -t134
              t189 = t8*t18*t25*zh*1.0e+1_r64
              t193 = -t141
              t194 = chi*t2*t8*t25*(4.0_r64*ic)
              t195 = chi*t3*t7*t26*(4.0_r64*ic)
              t203 = t82*zh*2.0_r64
              t204 = rhop*t2*t36*2.0_r64
              t206 = t2*t3*t8*t25*2.0_r64
              t207 = t2*t3*t7*t26*2.0_r64
              t216 = t2*t4*t26*(1.0_r64*ic)
              t217 = nr*t2**2*t26*(4.0_r64*ic)
              t218 = t3*t5*t25*(1.0_r64*ic)
              t219 = nrp*t3**2*t25*(4.0_r64*ic)
              t222 = -t160
              t223 = -t161
              t224 = rho*t7*t26*(8.0_r64*ic)
              t226 = -t169
              t227 = rhop*t8*t25*(8.0_r64*ic)
              t228 = -t175
              t229 = t2*t3*t8*t25*zh
              t230 = t3**2*t7*t26*zh
              t236 = -t185
              t237 = t2*t7*t26*(1.0_r64*ic)
              t238 = chi*t2*t7*t26*(4.0_r64*ic)
              t239 = t3*t8*t25*(1.0_r64*ic)
              t240 = chi*t3*t8*t25*(4.0_r64*ic)
              t241 = -t191
              t242 = t25*t26*zh*(8.0_r64*ic)
              t245 = t8+t26+t32
              t247 = nr*t2*t22*t26*3.0_r64
              t248 = nrp*t3*t22*t25*3.0_r64
              t249 = rhop*t3*t8*t25*zh*3.0_r64
              t251 = t3**2*t8*t25*zh*2.0_r64
              t258 = -t210
              t261 = t18*t36*zh*1.7e+1_r64
              t264 = -t213
              t269 = -t220
              t270 = -t221
              t272 = t2**2*t5*t25*zh*2.0_r64
              t273 = rhop*t2*t36*zh*3.0_r64
              t275 = rhop*t2**2*t8*t25*2.0_r64
              t277 = rhop*t2*t8*t25*zh*6.0_r64
              t282 = t3*t25*t26*(8.0_r64*ic)
              t283 = t2*t25*t26*(8.0_r64*ic)
              t284 = -t232
              t285 = -t233
              t288 = rhop*t2**2*t5*t25*1.2e+1_r64
              t289 = nz*rhop*t2**2*t5*t30
              t290 = nz*nzp*t2**2*t18*t30
              t294 = -t252
              t295 = -t254
              t296 = -t262
              t297 = -t263
              t298 = -t265
              t299 = -t266
              t301 = -t279
              t68 = t67**2
              t126 = t81*4.0_r64
              t143 = -t91
              t158 = -t108
              t165 = -t120
              t170 = -t122
              t171 = -t123
              t186 = -t133
              t196 = -t148
              t197 = -t149
              t202 = t81*zh*2.0_r64
              t225 = -t168
              t243 = -t194
              t244 = t4+t7+t38
              t255 = -t204
              t256 = -t207
              t267 = -t216
              t268 = -t218
              t271 = -t224
              t276 = t207*zh
              t286 = -t237
              t287 = -t240
              t293 = -t251
              t300 = -t277
              et1 = t45+t46+t47+t48-t57*1.24e+2_r64-t58*1.24e+2_r64+t69+t102+t126+t127+t180+t181+t183+t184+t236+t291+t292+t299+t18*t36*9.1e+1_r64+t2**2*t36*3.0_r64-t3**2*t37*3.0_r64+t2**2*t43*3.0_r64+t3**2*t43*3.0_r64-nr*t2**3*t26*2.0_r64-nrp*t3**3*t25*2.0_r64-rhop*t2*t36*6.5e+1_r64+t4*t5*t6*5.0e+1_r64-t4*t6*t8*3.3e+1_r64+t5*t6*t7*3.3e+1_r64-t6*t7*t8*2.1e+1_r64+rho*t2*t4*t8*2.9e+1_r64-rho*t2*t5*t7*9.1e+1_r64+rho*t2*t7*t8*2.9e+1_r64+rho*t2*t4*t26*1.9e+1_r64+rho*t2*t7*t26*1.9e+1_r64-rhop*t3*t5*t7*2.9e+1_r64+rhop*t3*t7*t8*2.9e+1_r64+rhop*t3*t5*t25*1.9e+1_r64-rhop*t3*t8*t25*1.9e+1_r64+t2*t3*t5*t7*6.5e+1_r64+t2*t3*t8*t25-t2*t3*t25*t26*5.2e+1_r64
              et2 = t3*t7*t26*t30-chi*t2**2*t7*t26*2.0_r64+chi*t2**2*t8*t25*2.0_r64-chi*t3**2*t7*t26*2.0_r64+chi*t3**2*t8*t25*2.0_r64-chi*t2**2*t25*t26*8.0_r64-chi*t3**2*t25*t26*8.0_r64+nr*rhop*t2**2*t5*1.78e+2_r64+t2*t3*t7*t8*t10-t2*t3*t7*t10*t26*4.0_r64+t2*t3*t8*t10*t25*4.0_r64+t2*t3*t10*t25*t26*1.6e+1_r64
              et3 = t45+t46+t47+t48-t57*1.26e+2_r64-t58*1.26e+2_r64+t69+t102+t126+t127+t180+t181+t183+t184+t236+t291+t292+t299+t18*t36*9.3e+1_r64+t2**2*t36*2.0_r64-t3**2*t37*2.0_r64+t2**2*t43*2.0_r64+t3**2*t43*2.0_r64-nr*t2**3*t26*3.0_r64-nrp*t3**3*t25*3.0_r64-rhop*t2*t36*6.6e+1_r64+t4*t5*t6*5.4e+1_r64-t4*t6*t8*3.4e+1_r64+t5*t6*t7*3.4e+1_r64-t6*t7*t8*2.2e+1_r64+rho*t2*t4*t8*3.0e+1_r64-rho*t2*t5*t7*9.3e+1_r64+rho*t2*t7*t8*3.0e+1_r64+rho*t2*t4*t26*2.0e+1_r64+rho*t2*t7*t26*2.0e+1_r64-rhop*t3*t5*t7*3.0e+1_r64+rhop*t3*t7*t8*3.0e+1_r64+rhop*t3*t5*t25*2.0e+1_r64-rhop*t3*t8*t25*2.0e+1_r64+t2*t3*t5*t7*6.6e+1_r64-t2*t3*t25*t26*4.6e+1_r64
              et4 = chi*t2**2*t7*t26*(-1.0_r64)*3.0_r64+chi*t2**2*t8*t25*4.0_r64-chi*t3**2*t7*t26*4.0_r64+chi*t3**2*t8*t25*3.0_r64-chi*t2**2*t25*t26*6.0_r64-chi*t3**2*t25*t26*6.0_r64+nr*rhop*t2**2*t5*1.72e+2_r64+t2*t3*t7*t8*t10*2.0_r64-t2*t3*t7*t10*t26*3.0_r64+t2*t3*t8*t10*t25*3.0_r64+t2*t3*t10*t25*t26*1.2e+1_r64
              et5 = -t109-t229*6.0_r64+chi*t230+t57*zh*9.1e+1_r64+t58*zh*2.9e+1_r64+t72*zh+t81*zh+t82*zh*3.0_r64+t3**2*t37*zh*3.0_r64-t3**2*t43*zh*3.0_r64-t18*t36*zh*5.8e+1_r64+nz*t5*t6*t15*3.0_r64+nz*t5*t6*t18*3.0_r64-nzp*t6*t7*t18*6.0_r64-nr*t5*t19*zh*1.5e+1_r64+rhop*t2*t36*zh*2.3e+1_r64-t4*t5*t6*zh*3.3e+1_r64+t4*t6*t8*zh*9.0_r64-t5*t6*t7*zh*2.1e+1_r64+t6*t7*t8*zh*9.0_r64+nr*nzp*t2**2*t19*8.0_r64-nz*t2**2*t5*t6*2.0_r64-nz*t2*t3**3*t26*2.0_r64+nz*t3**3*t5*t30+nzp*t2*t3**3*t7*2.0_r64+nrp*t3**3*t25*zh*2.0_r64-nr*rhop*t2**2*t5*zh*6.5e+1_r64-nz*rhop*t2*t5*t6*1.0e+1_r64+nzp*rhop*t2*t6*t7*8.0_r64
              et6 = nz*t2*t3*t5*t18*(-1.0_r64)*2.0_r64+nzp*t2*t3*t7*t18*4.0_r64+rho*t2*t5*t7*zh*2.9e+1_r64+rhop*t3*t5*t7*zh*2.9e+1_r64-rhop*t3*t7*t8*zh*2.9e+1_r64-rhop*t3*t5*t25*zh*1.9e+1_r64+rhop*t3*t8*t25*zh*1.9e+1_r64-t2*t3*t5*t7*zh*4.4e+1_r64+t2*t3*t7*t8*zh*2.3e+1_r64+t2*t3*t7*t26*zh*3.0_r64+t2*t3*t25*t26*zh*2.0_r64-chi*t3**2*t8*t25*zh*2.0_r64+chi*t3**2*t25*t26*zh*4.0_r64+nz*rhop*t2**2*t3*t5*1.2e+1_r64-nzp*rhop*t2**2*t3*t7*8.0_r64-nzp*rhop*t2**2*t3*t25*4.0_r64+nz*t2**2*t3*t5*t30
              et7 = t110+t230*6.0_r64+chi*t276-nzp*t4*t6**2*3.0_r64-nzp*t6**2*t7*6.0_r64-nzp*t4*t20*3.0_r64+nzp*t2*t4*t19*1.0e+1_r64+nzp*t2*t7*t19*8.0_r64+nzp*t2*t19*t25*2.0_r64-t4*t5*t15*zh*1.5e+1_r64-t4*t5*t18*zh*3.3e+1_r64-t5*t7*t15*zh*1.5e+1_r64+t4*t8*t18*zh*2.1e+1_r64-t5*t7*t18*zh*9.0_r64+t7*t8*t18*zh*9.0_r64-t2*t3*t43*zh*3.0_r64+t5*t18*t25*zh*2.9e+1_r64-t8*t18*t25*zh*2.9e+1_r64+nzp*rhop*t2**3*t25-nz*t2**2*t5*t18*8.0_r64+nzp*t2**2*t4*t18*2.0_r64+nzp*t2**2*t7*t18*4.0_r64-nzp*t2**2*t18*t25*1.2e+1_r64-t2**2*t4*t5*zh*1.9e+1_r64-t2**2*t5*t7*zh*1.9e+1_r64-t3**2*t5*t7*zh*2.3e+1_r64+t3**2*t7*t8*zh*2.3e+1_r64
              et8 = t2**2*t5*t25*zh+t3**2*t5*t25*zh*3.0_r64-t2**2*t8*t25*zh*3.0_r64-t3**2*t8*t25*zh*3.0_r64+t2**2*t25*t26*zh*2.0_r64+t3**2*t25*t26*zh*2.0_r64+nz*t2**2*t3**2*t5*4.0_r64+nzp*t2**2*t3**2*t7*2.0_r64+nzp*t2**2*t3**2*t25*2.0_r64+rhop*t2*t4*t5*zh*9.1e+1_r64-rhop*t2*t4*t8*zh*2.9e+1_r64+rhop*t2*t5*t7*zh*5.8e+1_r64-rhop*t2*t7*t8*zh*2.9e+1_r64-rhop*t2*t5*t25*zh*6.5e+1_r64+rhop*t2*t8*t25*zh*4.4e+1_r64-nzp*rhop*t2*t3**2*t7*8.0_r64+nzp*rhop*t2*t3**2*t25+chi*t2*t3*t25*t26*zh*4.0_r64+chi*t3*t8*t25*t30*zh
              et9 = t159-t231*1.6e+1_r64+t288+t3*t132+t3*t157+t3*t176+t3*t206+t3*t248+t4*t5*t19*3.0_r64-t4*t8*t19*6.0_r64-t3*t22*t37*2.3e+1_r64+t3*t22*t43*2.3e+1_r64+rho*t4*t5*t6*3.0_r64+rho*t5*t6*t7*6.0_r64-rhop*t4*t5*t22*2.1e+1_r64+rhop*t4*t8*t22*9.0_r64-rhop*t5*t7*t22*9.0_r64+rhop*t7*t8*t22*9.0_r64+rhop*t5*t22*t25*2.9e+1_r64-rhop*t8*t22*t25*2.9e+1_r64-t2*t4*t5*t18*1.0e+1_r64+t2*t4*t8*t18*8.0_r64-t2*t5*t7*t18*8.0_r64+t2*t4*t5*t22*2.9e+1_r64+t2*t7*t8*t18*1.6e+1_r64+t2*t5*t7*t22*2.9e+1_r64-t2*t5*t18*t25*2.0_r64+t2*t8*t18*t25*4.0_r64-t2*t5*t22*t25*4.4e+1_r64+t2*t8*t22*t25*2.3e+1_r64+t2*t22*t25*t26*3.0_r64+nz*nzp*t2**2*t19*8.0_r64-rhop*t2**2*t4*t5*2.0_r64
              et10 = rhop*t2**2*t5*t7*(-1.0_r64)*4.0_r64-rhop*t2**2*t8*t25*8.0_r64+t2*t3**2*t5*t7*8.0_r64-t2*t3**2*t7*t26*2.0_r64-t2*t3**2*t25*t26*2.0_r64+t2**2*t5*t25*t30-nz*nzp*rhop*t2**2*t3**2*4.0_r64-chi*t3*t8*t22*t25*3.0_r64
              pcsm(1,jp) = ipi*((t56*t68*(t61+t65+t71+t77+t87+t112+t119+t135+t139+t144+t145+t158+t159+t171+t178+t182+t187+t208+t231+t253+t257+t258+t275+t280+t290+t297+t298-t3*t81*8.0_r64-t2**3*t5*t25*8.0_r64-t4*t5*t19*2.8e+1_r64+t4*t8*t19*1.9e+1_r64-rho*t4*t5*t6*2.8e+1_r64-rho*t5*t6*t7*1.9e+1_r64+t2*t4*t5*t18*6.8e+1_r64+t2*t5*t18*t25*4.0_r64+t2*t8*t18*t25+nz*nzp*t2**2*t19*6.0_r64+rhop*t2**2*t4*t5*4.0_r64-rhop*t2**2*t5*t25*3.6e+1_r64-t2*t3**2*t7*t26*4.0_r64+t2*t3**2*t8*t25*4.0_r64+t2*t3**2*t25*t26*1.6e+1_r64-nz*nzp*rhop*t2**2*t3**2*2.0_r64+rhop*t2*t5*t7*t30))/8.0_r64)
              pcsm(2,jp) = ipi*(t56*t68*(t61+t65+t71+t77+t87+t112+t119+t135+t139+t144+t145+t158+t171+t178+t182+t187+t208+t231+t253+t257+t258+t275+t280+t288+t290+t297+t298+t3*t81*1.0e+1_r64+t2**3*t5*t25*1.0e+1_r64-t4*t5*t19*1.0e+1_r64+t4*t8*t19*1.0e+1_r64-rho*t4*t5*t6*1.0e+1_r64-rho*t5*t6*t7*1.0e+1_r64+t2*t4*t5*t18*4.4e+1_r64-t2*t5*t18*t25*3.2e+1_r64+t2*t8*t18*t25*1.9e+1_r64-nz*nzp*rhop*t6**2-rhop*t2**2*t4*t5*3.2e+1_r64-rhop*t2**2*t5*t7*1.9e+1_r64+t2*t3**2*t7*t26*5.0_r64-t2*t3**2*t8*t25*5.0_r64-t2*t3**2*t25*t26*8.0_r64+nz*nzp*rhop*t2**2*t3**2)*(-1.0_r64/2.0_r64))
              pcsm(3,jp) = ipi*((t28*t35*t42*t55*(et1+et2))/8.0_r64)
              pcsm(4,jp) = ipi*(t28*t35*t42*t55*(et3+et4)*(-1.0_r64/2.0_r64))
              pcsm(5,jp) = ipi*(t28*t33*t41*t55*t244*t245*(t6+rho*t30+rhop*t31+t2*t3)*2.0_r64)
              pcsm(6,jp) = ipi*(t55*t67*(rhop*t36*(-1.0_r64)*(2.0_r64*ic)+nr*t2**2*t5*(1.6e+1_r64*ic)-nr*t5*t18*(9.0_r64*ic)+nr*t8*t18*(9.0_r64*ic)-rho*t4*t5*(1.9e+1_r64*ic)+rho*t4*t8*(1.0_r64*ic)-rho*t5*t7*(1.0e+1_r64*ic)+rho*t7*t8*(1.0_r64*ic)+rho*t4*t26*(1.0_r64*ic)+rho*t7*t26*(1.0_r64*ic)+t3*t5*t7*(1.0_r64*ic)-t3*t7*t8*(1.0_r64*ic)+t3*t5*t25*(8.0_r64*ic)+t3*t7*t26*(8.0_r64*ic)-t3*t8*t25*(8.0_r64*ic)-t3*t25*t26*(2.0e+1_r64*ic)-nz*nzp*rhop*t6*(3.0_r64*ic)+nr*rhop*t2*t5*(2.3e+1_r64*ic)+nz*nzp*rhop*t2*t3*(3.0_r64*ic))*(-1.0_r64/4.0_r64))
              pcsm(7,jp) = ipi*(t39*t55*t67*t244*t245*(1.0_r64*ic))
              pcsm(8,jp) = ipi*(t33*t41*t55*(t62+t95+t98+t100+t154+t156+t167+t172+t195+t197+t219+t222+t227+t228+t234+t267+t270+t283+t285+t286+t287+rhop*t4*t5*(1.1e+1_r64*ic)+t2*t5*t25*(5.0_r64*ic)-chi*t3*t25*t26*(4.0_r64*ic))*(-1.0_r64/4.0_r64))
              pcsm(9,jp) = ipi*(t33*t41*t55*(t62+t95+t98+t100+t154+t156+t167+t172+t195+t197+t219+t222+t227+t228+t234+t267+t270+t283+t285+t286+t287+rhop*t4*t5*(8.0_r64*ic)+t2*t5*t25*(1.1e+1_r64*ic)-chi*t3*t25*t26*(7.0_r64*ic)))
              pcsm(10,jp) = ipi*((t56*t68*(t60+t74+t84-t110+t117+t131+t137+t138+t146+t147+t162+t165+t179+t189+t201+t215+t230+t241+t250+t259+t272+t274+t289+t293+t294+t295+t296+t300+t301+t4*t5*t18*zh*1.9e+1_r64+nz*t2**2*t5*t18*6.0_r64+t3**2*t25*t26*zh*4.0_r64-nz*t2**2*t3**2*t5*2.0_r64+rhop*t2*t5*t25*zh))/8.0_r64)
              pcsm(11,jp) = ipi*(t56*t68*(t60+t74+t84+t117+t131+t137+t138+t146+t147+t162+t165+t179+t189+t201+t215+t230+t241+t250+t259+t272+t274+t289+t293+t294+t295+t296+t300+t301+nz*t6**2*t32+t4*t5*t18*zh*1.0e+1_r64-t3**2*t25*t26*zh*5.0_r64+nz*t2**2*t3**2*t5+rhop*t2*t5*t25*zh*1.9e+1_r64)*(-1.0_r64/2.0_r64))
              pcsm(12,jp) = ipi*(t35*t42*t56*(et5+et6)*(-1.0_r64/8.0_r64))
              pcsm(13,jp) = ipi*(rhop*t33*t39*t41*t56*t244*t245*zh*(-1.0_r64/2.0_r64))
              pcsm(14,jp) = ipi*(t18*t56*t67*(nrp*t4*t15*(9.0_r64*ic)+nrp*t7*t15*(9.0_r64*ic)+rhop*t4*t5*(1.9e+1_r64*ic)-rhop*t4*t8*(1.0e+1_r64*ic)+rhop*t5*t7*(1.0_r64*ic)-rhop*t7*t8*(1.0_r64*ic)-rhop*t5*t25*(1.0_r64*ic)+rhop*t8*t25*(1.0_r64*ic)-t2*t4*t5*(2.3e+1_r64*ic)+t2*t4*t8*(1.0_r64*ic)-t2*t5*t7*(2.0_r64*ic)+t2*t7*t8*(1.0_r64*ic)-t2*t4*t26*(8.0_r64*ic)-t2*t5*t25*(1.6e+1_r64*ic)-t2*t7*t26*(8.0_r64*ic)+t2*t8*t25*(8.0_r64*ic)+t2*t25*t26*(2.0e+1_r64*ic)+nz*nzp*rho*t6*(3.0_r64*ic)-nz*nzp*rhop*t2**2*(3.0_r64*ic))*(-1.0_r64/4.0_r64))
              pcsm(15,jp) = ipi*(t18*t40*t56*t67*t244*t245*(-1.0_r64)*(1.0_r64*ic))
              pcsm(16,jp) = ipi*((t29*t33*t41*t54*(t63+t94+t107+t150+t152+t153+t173+t174+t196+t217+t223+t226+t235+t238+t239+t243+t268+t269+t271+t282+t284+nr*t2**2*t5*(5.0_r64*ic)+rho*t4*t5*(1.1e+1_r64*ic)-chi*t2*t25*t26*(4.0_r64*ic)))/4.0_r64)
              pcsm(17,jp) = ipi*(-t29*t33*t41*t54*(t63+t94+t107+t150+t152+t153+t173+t174+t196+t217+t223+t226+t235+t238+t239+t243+t268+t269+t271+t282+t284+nr*t2**2*t5*(1.1e+1_r64*ic)+rho*t4*t5*(8.0_r64*ic)-chi*t2*t25*t26*(7.0_r64*ic)))
              pcsm(18,jp) = ipi*((t28*t54*t67*(t36+t43+t53-chi*t7*t26*4.0_r64+chi*t8*t25*4.0_r64+chi*t25*t26*1.6e+1_r64+nz*nzp*t6*2.0_r64+nr*rhop*t5*5.0_r64-nr*rhop*t8*5.0_r64+nrp*rho*t4*5.0_r64+nrp*rho*t7*5.0_r64-nr*t2*t5*1.8e+1_r64-nr*t2*t26*4.0_r64-nrp*t3*t25*4.0_r64-nz*nzp*t2*t3*2.0_r64))/8.0_r64)
              pcsm(19,jp) = ipi*(t28*t54*t67*(t36+t43+t53+chi*t7*t26*5.0_r64-chi*t8*t25*5.0_r64-chi*t25*t26*8.0_r64-nz*nzp*t6-nr*rhop*t5*4.0_r64+nr*rhop*t8*4.0_r64-nrp*rho*t4*4.0_r64-nrp*rho*t7*4.0_r64+nr*t2*t5*6.0_r64+nr*t2*t26*5.0_r64+nrp*t3*t25*5.0_r64+nz*nzp*t2*t3)*(-1.0_r64/2.0_r64))
              pcsm(20,jp) = ipi*((rhop*t33*t41*t55*(t44+t113+chi*t36+chi*t53+t4*t5*6.0_r64+t4*t8*3.0_r64-t5*t7*3.0_r64+t7*t8*3.0_r64+t4*t26*8.0_r64+t5*t25*8.0_r64+t7*t26*8.0_r64-t8*t25*8.0_r64-t25*t26*3.0e+1_r64+t7*t8*t10-t7*t10*t26*4.0_r64+t8*t10*t25*4.0_r64+t10*t25*t26*1.6e+1_r64-chi*nr*t2*t26*4.0_r64-chi*nrp*t3*t25*4.0_r64))/8.0_r64)
              pcsm(21,jp) = ipi*(rhop*t33*t41*t55*(t44+t113+chi*t36*2.0_r64-chi*t37*2.0_r64+t4*t5*4.0_r64+t4*t8*2.0_r64-t5*t7*2.0_r64+t7*t8*2.0_r64+t4*t26*7.0_r64+t5*t25*7.0_r64+t7*t26*7.0_r64-t8*t25*7.0_r64-t25*t26*2.4e+1_r64+t7*t8*t10*2.0_r64-t7*t10*t26*3.0_r64+t8*t10*t25*3.0_r64+t10*t25*t26*1.2e+1_r64-chi*nr*t2*t26*3.0_r64-chi*nrp*t3*t25*3.0_r64)*(-1.0_r64/2.0_r64))
              pcsm(22,jp) = ipi*((rhop*t55*t244*t245*2.0_r64)/t27)
              pcsm(23,jp) = ipi*(t29*t54*t67*(t64+t96+t151+t155+t163+t164+t166+t225+t242-nz*t2**2*t5*(3.0_r64*ic)+nz*t5*t15*(3.0_r64*ic))*(-1.0_r64/4.0_r64))
              pcsm(24,jp) = ipi*(t29*t54*t67*t244*t245*zh*(1.0_r64*ic))
              pcsm(25,jp) = ipi*((t28*t33*t41*t55*(t3*t168+t63*zh+t173*zh+rhop*t36*zh*(4.0_r64*ic)-nz*rhop*t5*t6*(1.0_r64*ic)+nzp*rhop*t6*t7*(2.0_r64*ic)+nz*t2*t5*t6*(1.0_r64*ic)-nr*t8*t18*zh*(5.0_r64*ic)+rho*t4*t5*zh*(4.0_r64*ic)+rho*t5*t7*zh*(4.0_r64*ic)-t3*t5*t7*zh*(4.0_r64*ic)-t3*t5*t25*zh*(1.0_r64*ic)+t3*t25*t26*zh*(4.0_r64*ic)-nz*t2**2*t3*t5*(1.0_r64*ic)+nz*rhop*t2*t3*t5*(1.0_r64*ic)-nzp*rhop*t2*t3*t7*(2.0_r64*ic)-nr*rhop*t2*t5*zh*(1.2e+1_r64*ic)))/4.0_r64)
              pcsm(26,jp) = ipi*(t56*t68*(t49+t78-t80+t85+t109-t116+t130-t136-t140-t142-t188+t190+t192+t198+t202+t203-t211-t212-t214-t249-t260+t261-t273+t276+t278+t281+nr*nzp*t20*t30+t4*t5*t6*zh*1.9e+1_r64-nr*nzp*t2**2*t19*6.0_r64+nr*rhop*t2**2*t5*zh+nzp*rhop*t4*t6*t30+t2*t3*t25*t26*zh*4.0_r64+t3*t8*t25*t30*zh+nzp*rhop*t2**2*t3*t25*2.0_r64)*(-1.0_r64/8.0_r64))
              pcsm(27,jp) = ipi*(t56*t68*(-t49+t59-t78+t80+t83+t116-t130+t136+t140+t142+t188-t190-t192-t202-t203+t211+t212+t214+t229+t249+t260-t261+t273-t278-t281-nzp*t4*t6*t18-t4*t5*t6*zh*1.0e+1_r64-nr*rhop*t2**2*t5*zh*1.9e+1_r64+nzp*t3*t18*t25*t30-t2*t3*t7*t26*zh*2.0_r64+t2*t3*t25*t26*zh*5.0_r64+nr*nzp*t2**2*t18*t30+nzp*rhop*t2**2*t3*t25)*(-1.0_r64/2.0_r64))
              pcsm(28,jp) = ipi*((t35*t42*t56*(et7+et8))/8.0_r64)
              pcsm(29,jp) = ipi*((rhop*t33*t40*t41*t56*t244*t245*zh)/2.0_r64)
              pcsm(30,jp) = ipi*(t55*t67*(t64+t96+t151+t155+t163+t164+t166+t225+t242-nzp*t4*t18*(3.0_r64*ic)+nzp*rhop*t2*t25*(3.0_r64*ic))*(-1.0_r64/4.0_r64))
              pcsm(31,jp) = ipi*(t55*t67*t244*t245*zh*(1.0_r64*ic))
              pcsm(32,jp) = ipi*((t33*t41*t55*(t36*zh*(4.0_r64*ic)-t37*zh*(4.0_r64*ic)+t43*zh*(4.0_r64*ic)+nzp*t4*t6*(1.0_r64*ic)+nzp*t6*t7*(2.0_r64*ic)-nr*nzp*t2*t18*(1.0_r64*ic)-chi*t7*t26*zh*(1.0_r64*ic)+chi*t25*t26*zh*(4.0_r64*ic)-nzp*t2*t3*t7*(2.0_r64*ic)+nzp*t2*t3*t25*(1.0_r64*ic)+nr*rhop*t5*zh*(4.0_r64*ic)-nr*rhop*t8*zh*(4.0_r64*ic)+nrp*rho*t4*zh*(5.0_r64*ic)+nrp*rho*t7*zh*(5.0_r64*ic)-nr*t2*t5*zh*(1.2e+1_r64*ic)-nr*t2*t26*zh*(1.0_r64*ic)-nr*nzp*rhop*t2**2*(1.0_r64*ic)))/4.0_r64)
              pcsm(33,jp) = ipi*((t28*t55*t68*(t50+t57+t58+t66+t70+t75+t79+t92+t114+t115+t121+t124+t132+t143+t157+t170+t176+t177+t186+t193+t205+t206+t247+t248+t255+t256+t264-t4*t5*t6*4.0_r64-t2*t3*t25*t26*2.0_r64+nr*rhop*t2**2*t5*6.0_r64))/8.0_r64)
              pcsm(34,jp) = ipi*(t28*t55*t68*(t50+t57+t58+t66+t70+t75+t79+t92+t114+t115+t121+t124+t132+t143+t157+t170+t176+t177+t186+t193+t205+t206+t247+t248+t255+t256+t264+t4*t6*t32+t2*t3*t25*t26)*(-1.0_r64/2.0_r64))
              pcsm(35,jp) = ipi*(t35*t42*t56*(et9+et10)*(-1.0_r64/8.0_r64))
              pcsm(36,jp) = ipi*(rhop*t22*t33*t41*t56*t244*t245*(-1.0_r64/2.0_r64))
            end block
          end do
          do md = 0, M
            rn = real(md,r64); n2 = rn*rn; ifn2 = 1.0_r64/(4.0_r64*n2-1.0_r64); n4 = n2*n2
            do jp = 1, p
              vk = vkmat(jp,md); ve = vemat(jp,md)
              coK(1)=pcsm(1,jp)+n2*pcsm(2,jp);   coK(3)=pcsm(10,jp)+n2*pcsm(11,jp); coK(5)=pcsm(18,jp)+n2*pcsm(19,jp)
              coK(7)=pcsm(26,jp)+n2*pcsm(27,jp); coK(9)=pcsm(33,jp)+n2*pcsm(34,jp)
              coK(2)=rn*(pcsm(6,jp)+n2*pcsm(7,jp));   coK(4)=rn*(pcsm(14,jp)+n2*pcsm(15,jp))
              coK(6)=rn*(pcsm(23,jp)+n2*pcsm(24,jp)); coK(8)=rn*(pcsm(30,jp)+n2*pcsm(31,jp))
              coE(1)=ifn2*(pcsm(3,jp)+n2*pcsm(4,jp)+n4*pcsm(5,jp)); coE(5)=ifn2*(pcsm(20,jp)+n2*pcsm(21,jp)+n4*pcsm(22,jp))
              coE(3)=pcsm(12,jp)+n2*pcsm(13,jp); coE(7)=pcsm(28,jp)+n2*pcsm(29,jp); coE(9)=pcsm(35,jp)+n2*pcsm(36,jp)
              coE(2)=rn*ifn2*(pcsm(8,jp)+n2*pcsm(9,jp)); coE(4)=rn*ifn2*(pcsm(16,jp)+n2*pcsm(17,jp))
              coE(6)=rn*pcsm(25,jp); coE(8)=rn*pcsm(32,jp)
              Be_rr(jp) = (coK(1)*vk+coE(1)*ve)*wsp(jp);  Be_rt(jp) = (coK(2)*vk+coE(2)*ve)*wsp(jp);  Be_rz(jp) = (coK(3)*vk+coE(3)*ve)*wsp(jp)
              Be_tr(jp) = (coK(4)*vk+coE(4)*ve)*wsp(jp);  Be_tt(jp) = (coK(5)*vk+coE(5)*ve)*wsp(jp);  Be_tz(jp) = (coK(6)*vk+coE(6)*ve)*wsp(jp)
              Be_zr(jp) = (coK(7)*vk+coE(7)*ve)*wsp(jp);  Be_zt(jp) = (coK(8)*vk+coE(8)*ve)*wsp(jp);  Be_zz(jp) = (coK(9)*vk+coE(9)*ve)*wsp(jp)
            end do
            do l = 1, p                                        ! fold each named block (sum over subpanel nodes via Lc) into A
              A(0*nt+nidx(i), 0*nso+cols+l, md+1) = A(0*nt+nidx(i), 0*nso+cols+l, md+1) + sum(Be_rr(1:p)*Lc(:,l))
              A(0*nt+nidx(i), 1*nso+cols+l, md+1) = A(0*nt+nidx(i), 1*nso+cols+l, md+1) + sum(Be_rt(1:p)*Lc(:,l))
              A(0*nt+nidx(i), 2*nso+cols+l, md+1) = A(0*nt+nidx(i), 2*nso+cols+l, md+1) + sum(Be_rz(1:p)*Lc(:,l))
              A(1*nt+nidx(i), 0*nso+cols+l, md+1) = A(1*nt+nidx(i), 0*nso+cols+l, md+1) + sum(Be_tr(1:p)*Lc(:,l))
              A(1*nt+nidx(i), 1*nso+cols+l, md+1) = A(1*nt+nidx(i), 1*nso+cols+l, md+1) + sum(Be_tt(1:p)*Lc(:,l))
              A(1*nt+nidx(i), 2*nso+cols+l, md+1) = A(1*nt+nidx(i), 2*nso+cols+l, md+1) + sum(Be_tz(1:p)*Lc(:,l))
              A(2*nt+nidx(i), 0*nso+cols+l, md+1) = A(2*nt+nidx(i), 0*nso+cols+l, md+1) + sum(Be_zr(1:p)*Lc(:,l))
              A(2*nt+nidx(i), 1*nso+cols+l, md+1) = A(2*nt+nidx(i), 1*nso+cols+l, md+1) + sum(Be_zt(1:p)*Lc(:,l))
              A(2*nt+nidx(i), 2*nso+cols+l, md+1) = A(2*nt+nidx(i), 2*nso+cols+l, md+1) + sum(Be_zz(1:p)*Lc(:,l))
            end do
          end do
        end do
      end do
    end do
    deallocate(tin, joa, ci, zc, zcn, nidx, znear, znearn, ikq, ikc, C1, C2, C3, C4, C5, As, Ad, A1, A2, A3, A4, Gcq, Gpc, ff)
  end subroutine axissymstok_dlpn_blockmat_nmode_r64
end module axissymstok_specialquad_mod
