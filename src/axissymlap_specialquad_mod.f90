module axissymlap_specialquad_mod
  use axilaplace3d_mod, only: r64, gauss_r64, lagrange_interp_r64
  use axisym_modal_green_mod, only: modal_green_r64, modal_green_all_r64, modal_green_all_far_r64
  use specialquad_mod, only: sdspecialquad_r64
  use axissymlap_kernelsplit_mod, only: axissymlap_slp_coef_r64, axissymlap_slp_coef_nmode_r64, &
       axissymlap_slpn_coef_r64, axissymlap_slpn_coef_nmode_r64, &
       axissymlap_dlp_coef_r64, axissymlap_dlp_coef_nmode_r64, &
       axissymlap_dlpn_coef_r64, axissymlap_dlpn_coef_nmode_r64
  implicit none
  private
  public :: axissymlap_slp_blockmat_r64
  public :: axissymlap_slp_blockmat_nmode_r64
  public :: axissymlap_slpn_blockmat_r64
  public :: axissymlap_slpn_blockmat_nmode_r64
  public :: axissymlap_dlp_blockmat_r64
  public :: axissymlap_dlp_blockmat_nmode_r64
  public :: axissymlap_dlpn_blockmat_r64
  public :: axissymlap_dlpn_blockmat_nmode_r64
contains

  subroutine axissymlap_slp_blockmat_r64(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, A)
    ! scalar Laplace SLP block (mode 0), real A(nt, np*p) -- mirror of LapSLPAxiBlockMat_nmode at n=0.
    integer(8),   intent(in)    :: nt, p, np, iside, iclosed
    complex(r64), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: sws(p*np), tpan(np+1)
    real(r64),    intent(inout) :: A(nt, np*p)
    integer(8) :: q, i, j, jp, iq, ia, l, nc, cjo, k, jo, ne, npa
    real(r64)  :: twopi, sumws, rho, rhop, zh, rr2, chi, rt, zti, vk, ve, Fn, An, dFn, SK
    real(r64)  :: t1, tN, tm, denom, rlo, rhi, tgi, split, ff(p)
    complex(r64) :: ic, za, zb
    real(r64)    :: tglp(p), wglp(p), Dp(p,p), tglq(2*p), wglq(2*p), Dq(2*p,2*p), IP2(2*p,p)
    real(r64)    :: rk(p), re2(2), IPk(p,p), IPe2(2,p), wsq(2*p), wsp(p)
    complex(r64) :: Yp(p), Yq(2*p), dYp(p), dYq(2*p), nvq(2*p), wxpq(2*p), ec2(2), xc(p,np)
    real(r64),    allocatable :: tin(:)
    logical,      allocatable :: cl(:)
    integer(8),   allocatable :: ci(:)
    complex(r64), allocatable :: zc(:), Ad(:,:)
    real(r64),    allocatable :: C1(:,:), C2(:,:), As(:,:), A1(:,:), A2(:,:), A3(:,:), A4(:,:), Gcq(:,:), Gpc(:,:)
    q = 2*p; ic = (0.0_r64,1.0_r64); twopi = 2.0_r64*acos(-1.0_r64)
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
    ! dyadic pole refinement of tpan -> tin (mirror LapSLPAxiBlockMat_nmode, refine 1st/last)
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
    allocate(cl(nt), ci(nt), zc(nt), C1(nt,q), C2(nt,q))
    allocate(As(q,nt), Ad(q,nt), A1(q,nt), A2(q,nt), A3(q,nt), A4(q,nt), Gcq(nt,q), Gpc(nt,p))
    A = 0.0_r64
    do k = 1, npa
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
      Yq = matmul(IP2, Yp); dYq = matmul(Dq, Yq)
      do iq = 1, q
        nvq(iq) = -ic*dYq(iq)/abs(dYq(iq)); wsq(iq) = wglq(iq)*abs(dYq(iq)); wxpq(iq) = dYq(iq)*wglq(iq)
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
      ! near: 0th SLP close (coef + sdspecialquad log), fold q->p (IP2) -> coarse (IPk)
      if (nc > 0) then
        call axissymlap_slp_coef_r64(nc, zc(1:nc), q, Yq, C1(1:nc,1:q), C2(1:nc,1:q))
        call sdspecialquad_r64(nc, zc(1:nc), q, Yq, nvq, wxpq, za, zb, iside, &
             As(1:q,1:nc), Ad(1:q,1:nc), A1(1:q,1:nc), A2(1:q,1:nc), A3(1:q,1:nc), A4(1:q,1:nc))
        do iq = 1, q
          do ia = 1, nc
            Gcq(ia,iq) = twopi*C1(ia,iq)*As(iq,ia) + C2(ia,iq)*wsq(iq)
          end do
        end do
        Gpc(1:nc,1:p) = matmul(Gcq(1:nc,1:q), IP2)
        do ia = 1, nc
          do l = 1, p
            A(ci(ia), cjo+l) = A(ci(ia), cjo+l) + sum(Gpc(ia,1:p)*IPk(:,l))
          end do
        end do
      end if
      ! far: analytic 0th single carrier K_0 = SK*VK on the refined p nodes, fold via IPk
      do i = 1, nt
        if (cl(i)) cycle
        rt = real(tx(i),r64); zti = aimag(tx(i))
        do jp = 1, p
          rhop = real(Yp(jp),r64); zh = zti - aimag(Yp(jp)); rho = rt
          rr2 = (rho-rhop)**2 + zh*zh; chi = 1.0_r64 + rr2/(2.0_r64*rho*rhop)
          call modal_green_r64(chi, 0_8, vk, ve, Fn, An, dFn)
          SK = sqrt(rhop/rho)/twopi
          ff(jp) = SK*vk*wsp(jp)
        end do
        do l = 1, p
          A(i, cjo+l) = A(i, cjo+l) + sum(ff(1:p)*IPk(:,l))
        end do
      end do
    end do
    deallocate(tin, cl, ci, zc, C1, C2, As, Ad, A1, A2, A3, A4, Gcq, Gpc)
  end subroutine axissymlap_slp_blockmat_r64

  subroutine axissymlap_slp_blockmat_nmode_r64(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, A)
    integer(8),   intent(in)    :: nt, p, np, M, iside, iclosed
    complex(r64), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: sws(p*np), tpan(np+1)
    real(r64),    intent(inout) :: A(nt, np*p, M+1)
    integer(8) :: q, nsa, i, j, jp, iq, ia, l, nc, cjo, md, k, jo, ne, npa, coff, K1, K2, nr1, nr2, c1b, c2b
    real(r64)  :: twopi, sumws, rho, rhop, zh, rr2, chi, rt, zti
    real(r64)  :: t1, tN, tm, denom, rlo, rhi, tgi, split, ff(p)
    complex(r64) :: ic, za, zb
    real(r64)    :: tglp(p), wglp(p), Dp(p,p), tglq(2*p), wglq(2*p), Dq(2*p,2*p), IP2(2*p,p)
    real(r64)    :: rk(p), re2(2), IPk(p,p), IPe2(2,p), wsq(2*p), wsp(p)
    complex(r64) :: Yp(p), Yq(2*p), dYp(p), dYq(2*p), nvq(2*p), wxpq(2*p), ec2(2), xc(p,np)
    real(r64),    allocatable :: tin(:), L1(:,:), L2(:,:), t_aux1(:), t_aux2(:), Be(:,:,:)
    logical,      allocatable :: cl(:)
    integer(8),   allocatable :: ci(:)
    complex(r64), allocatable :: zc(:), Ad(:,:)
    real(r64),    allocatable :: C1(:,:,:), C2(:,:,:), As(:,:), A1(:,:), A2(:,:), A3(:,:), A4(:,:), Gcq(:,:), Gpc(:,:)
    real(r64)    :: vka(0:M), vea(0:M), vkmat(p,0:M), SKw(p)
    q = 2*p; ic = (0.0_r64,1.0_r64); twopi = 2.0_r64*acos(-1.0_r64)
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
    nr1 = 0
    do while (tin(2) > t1)
      split = 0.5_r64*(tin(1)+tin(2))
      do i = ne, 2, -1; tin(i+1) = tin(i); end do
      tin(2) = split; ne = ne+1; nr1 = nr1+1
    end do
    nr2 = 0
    do while (tin(ne-1) < tN)
      split = 0.5_r64*(tin(ne)+tin(ne-1))
      tin(ne+1) = tin(ne); tin(ne) = split; ne = ne+1; nr2 = nr2+1
    end do
    npa = ne-1; nsa = npa*p
    K1 = nr1+1; K2 = nr2+1; c1b = K1*p; c2b = (K1+np-2)*p   ! Be col blocks: aux1=1:c1b, mid=c1b+1:c2b, aux2=c2b+1:nsa
    allocate(cl(nt), ci(nt), zc(nt), C1(nt,q,M+1), C2(nt,q,M+1))
    allocate(As(q,nt), Ad(q,nt), A1(q,nt), A2(q,nt), A3(q,nt), A4(q,nt), Gcq(nt,q), Gpc(nt,p))
    allocate(L1(K1*p,p), L2(K2*p,p), t_aux1(K1*p), t_aux2(K2*p), Be(nt,nsa,M+1))
    do j = 1, K1                                          ! pole-1 sub-panels (parent panel 1), nodes in parent-1 coords
      do i = 1, p
        tgi = tin(j) + (1.0_r64+tglp(i))/2.0_r64*(tin(j+1)-tin(j))
        t_aux1((j-1)*p+i) = (2.0_r64*tgi - (tpan(1)+tpan(2)))/(tpan(2)-tpan(1))
      end do
    end do
    do j = 1, K2                                          ! pole-2 sub-panels (parent panel np), nodes in parent-np coords
      k = (K1+np-2) + j
      do i = 1, p
        tgi = tin(k) + (1.0_r64+tglp(i))/2.0_r64*(tin(k+1)-tin(k))
        t_aux2((j-1)*p+i) = (2.0_r64*tgi - (tpan(np)+tpan(np+1)))/(tpan(np+1)-tpan(np))
      end do
    end do
    call lagrange_interp_r64(p, tglp, K1*p, t_aux1, L1)   ! == MATLAB  L1 = interpmat_1d(t_aux1, gauss(p))
    call lagrange_interp_r64(p, tglp, K2*p, t_aux2, L2)
    Be = 0.0_r64
    do k = 1, npa
      tm = 0.5_r64*(tin(k)+tin(k+1)); jo = 1
      do j = 1, np
        if (tm >= tpan(j) .and. tm < tpan(j+1)) jo = j
      end do
      coff = (k-1)*p
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
      Yq = matmul(IP2, Yp); dYq = matmul(Dq, Yq)
      do iq = 1, q
        nvq(iq) = -ic*dYq(iq)/abs(dYq(iq)); wsq(iq) = wglq(iq)*abs(dYq(iq)); wxpq(iq) = dYq(iq)*wglq(iq)
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
      if (nc > 0) then
        call axissymlap_slp_coef_nmode_r64(nc, zc(1:nc), q, Yq, M, C1(1:nc,1:q,:), C2(1:nc,1:q,:))
        call sdspecialquad_r64(nc, zc(1:nc), q, Yq, nvq, wxpq, za, zb, iside, &
             As(1:q,1:nc), Ad(1:q,1:nc), A1(1:q,1:nc), A2(1:q,1:nc), A3(1:q,1:nc), A4(1:q,1:nc))
        do md = 0, M
          do iq = 1, q
            do ia = 1, nc
              Gcq(ia,iq) = twopi*C1(ia,iq,md+1)*As(iq,ia) + C2(ia,iq,md+1)*wsq(iq)
            end do
          end do
          Gpc(1:nc,1:p) = matmul(Gcq(1:nc,1:q), IP2)         ! 2p -> p (refined nodes)
          do ia = 1, nc                                      ! write into Be (refined nodes, NO IPk)
            do l = 1, p
              Be(ci(ia), coff+l, md+1) = Be(ci(ia), coff+l, md+1) + Gpc(ia,l)
            end do
          end do
        end do
      end if
      do i = 1, nt
        if (cl(i)) cycle
        rt = real(tx(i),r64); zti = aimag(tx(i))
        do jp = 1, p
          rhop = real(Yp(jp),r64); zh = zti - aimag(Yp(jp)); rho = rt
          rr2 = (rho-rhop)**2 + zh*zh; chi = 1.0_r64 + rr2/(2.0_r64*rho*rhop)
          call modal_green_all_far_r64(chi, M, vka, vea)
          vkmat(jp,0:M) = vka(0:M); SKw(jp) = sqrt(rhop/rho)/twopi*wsp(jp)
        end do
        do md = 0, M
          do jp = 1, p
            ff(jp) = SKw(jp)*vkmat(jp,md)
          end do
          do jp = 1, p                                       ! write into Be at refined col coff+jp (NO IPk)
            Be(i, coff+jp, md+1) = Be(i, coff+jp, md+1) + ff(jp)
          end do
        end do
      end do
    end do
    ! ===== refined -> coarse projection:  A = Be * intMat,  intMat = blkdiag(L1, eye, L2) =====
    A = 0.0_r64
    do md = 0, M
      A(1:nt, 1:p,             md+1) = A(1:nt, 1:p,             md+1) + matmul(Be(:,1:c1b,md+1),   L1)
      A(1:nt, p+1:(np-1)*p,    md+1) = A(1:nt, p+1:(np-1)*p,    md+1) +        Be(:,c1b+1:c2b,md+1)
      A(1:nt, (np-1)*p+1:np*p, md+1) = A(1:nt, (np-1)*p+1:np*p, md+1) + matmul(Be(:,c2b+1:nsa,md+1), L2)
    end do
    deallocate(tin, L1, L2, t_aux1, t_aux2, Be, cl, ci, zc, C1, C2, As, Ad, A1, A2, A3, A4, Gcq, Gpc)
  end subroutine axissymlap_slp_blockmat_nmode_r64

  subroutine axissymlap_slpn_blockmat_r64(nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, A)
    ! scalar Laplace SLPn (S'=d_n S, target-normal traction) block (mode 0), real A(nt, np*p).
    integer(8),   intent(in)    :: nt, p, np, iside, iclosed
    complex(r64), intent(in)    :: tx(nt), tnx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: sws(p*np), tpan(np+1)
    real(r64),    intent(inout) :: A(nt, np*p)
    integer(8) :: q, i, j, jp, ia, l, nc, cjo, k, jo, ne, npa
    real(r64)  :: twopi, sumws, rho, rhop, zh, rr2, chi, rt, zti, vk, ve, Fn, An, dFn, SK, nr, nz, Qd, rn, cv
    real(r64)  :: tfirst, tlast, t0p, t1p, gam, denom, rlo, rhi, tgi, split, ff(p)
    logical    :: changed
    complex(r64) :: ic, za, zb
    real(r64)    :: tglp(p), wglp(p), Dp(p,p), IPk(p,p), IPe2(2,p), rk(p), re2(2), wsp(p)
    complex(r64) :: Yp(p), dYp(p), nvp(p), wxpp(p), ec2(2), xc(p,np)
    real(r64),    allocatable :: tin(:)
    logical,      allocatable :: cl(:)
    integer(8),   allocatable :: ci(:)
    complex(r64), allocatable :: zc(:), zcn(:), Ad(:,:)
    real(r64),    allocatable :: C1(:,:), C2(:,:), C3(:,:), As(:,:), A1(:,:), A2(:,:), A3(:,:), A4(:,:), Gcq(:,:)
    q = 2*p; ic = (0.0_r64,1.0_r64); twopi = 2.0_r64*acos(-1.0_r64)
    call gauss_r64(p, tglp, wglp, Dp)
    do j = 1, np
      do i = 1, p
        xc(i,j) = sx((j-1)*p+i)
      end do
    end do
    ! ---- graded dyadic pole refinement of tpan -> tin (gam=0.4 on every panel) ----
    tfirst = tpan(1)  + (1.0_r64+tglp(1))/2.0_r64*(tpan(2)-tpan(1))
    tlast  = tpan(np) + (1.0_r64+tglp(p))/2.0_r64*(tpan(np+1)-tpan(np))
    t0p = tpan(1); t1p = tpan(np+1)
    allocate(tin(np+1+1024)); ne = np+1; tin(1:ne) = tpan(1:np+1)
    do while (tin(2) > tfirst)
      split = 0.5_r64*(tin(1)+tin(2))
      do i = ne, 2, -1; tin(i+1) = tin(i); end do
      tin(2) = split; ne = ne+1
    end do
    do while (tin(ne-1) < tlast)
      split = 0.5_r64*(tin(ne)+tin(ne-1))
      tin(ne+1) = tin(ne); tin(ne) = split; ne = ne+1
    end do
    gam = 0.4_r64; changed = .true.
    do while (changed)
      changed = .false.; i = 1
      do while (i < ne)
        denom = min(tin(i)-t0p, t1p-tin(i+1))
        if (denom > 0.0_r64 .and. (tin(i+1)-tin(i)) > gam*denom) then
          do l = ne, i+1, -1; tin(l+1) = tin(l); end do
          tin(i+1) = 0.5_r64*(tin(i)+tin(i+1)); ne = ne+1; changed = .true.; i = i+2
        else
          i = i+1
        end if
      end do
    end do
    npa = ne-1
    allocate(cl(nt), ci(nt), zc(nt), zcn(nt), C1(nt,q), C2(nt,q), C3(nt,q))
    allocate(As(q,nt), Ad(q,nt), A1(q,nt), A2(q,nt), A3(q,nt), A4(q,nt), Gcq(nt,q))
    A = 0.0_r64
    do k = 1, npa
      tgi = 0.5_r64*(tin(k)+tin(k+1)); jo = 1
      do j = 1, np
        if (tgi >= tpan(j) .and. tgi < tpan(j+1)) jo = j
      end do
      denom = tpan(jo+1) - tpan(jo)
      do i = 1, p
        rlo = tin(k) + (1.0_r64+tglp(i))/2.0_r64*(tin(k+1)-tin(k))
        rk(i) = (2.0_r64*rlo - (tpan(jo)+tpan(jo+1)))/denom
      end do
      call lagrange_interp_r64(p, tglp, p, rk, IPk)
      Yp = matmul(IPk, xc(:,jo))
      rlo = (2.0_r64*tin(k)   - (tpan(jo)+tpan(jo+1)))/denom
      rhi = (2.0_r64*tin(k+1) - (tpan(jo)+tpan(jo+1)))/denom
      re2(1) = rlo; re2(2) = rhi
      call lagrange_interp_r64(p, tglp, 2_8, re2, IPe2)
      ec2 = matmul(IPe2, xc(:,jo)); za = ec2(1); zb = ec2(2)
      cjo = (jo-1)*p
      dYp = matmul(Dp, Yp)
      do jp = 1, p
        nvp(jp) = -ic*dYp(jp)/abs(dYp(jp)); wsp(jp) = wglp(jp)*abs(dYp(jp)); wxpp(jp) = dYp(jp)*wglp(jp)
      end do
      sumws = sum(wsp)
      nc = 0
      do i = 1, nt
        cl(i) = (abs(tx(i)-za) + abs(tx(i)-zb)) < 1.5_r64*sumws
        if (cl(i)) then; nc = nc + 1; ci(nc) = i; zc(nc) = tx(i); zcn(nc) = tnx(i); end if
      end do
      ! near: p-DIRECT 0th S' close on the aux nodes (coef + sdspecialquad As, target Cauchy -Ad/nvp)
      ! + analytic self-node override at the coincident pair, fold to coarse via IPk (== reference)
      if (nc > 0) then
        call axissymlap_slpn_coef_r64(nc, zc(1:nc), zcn(1:nc), p, Yp, C1(1:nc,1:p), C2(1:nc,1:p), C3(1:nc,1:p))
        call sdspecialquad_r64(nc, zc(1:nc), p, Yp, nvp, wxpp, za, zb, iside, &
             As(1:p,1:nc), Ad(1:p,1:nc), A1(1:p,1:nc), A2(1:p,1:nc), A3(1:p,1:nc), A4(1:p,1:nc))
        rn = 2.5_r64*log(2.0_r64)
        do jp = 1, p
          do ia = 1, nc
            cv = C2(ia,jp)
            if (abs(zc(ia)-Yp(jp)) < 1.0e-13_r64) then                   ! coincident self node: analytic c2 limit
              rho = real(zc(ia),r64); nr = real(zcn(ia),r64)
              cv = nr/(2.0_r64*twopi*rho)*(1.0_r64 - rn - 0.5_r64*log(2.0_r64*rho*rho))
            end if
            Gcq(ia,jp) = twopi*C1(ia,jp)*As(jp,ia) + cv*wsp(jp) &
                       - twopi*C3(ia,jp)*real(zcn(ia)*Ad(jp,ia)/nvp(jp), r64)   ! Sgrad = -Ad/nvp
          end do
        end do
        do ia = 1, nc
          do l = 1, p
            A(ci(ia), cjo+l) = A(ci(ia), cjo+l) + sum(Gcq(ia,1:p)*IPk(:,l))
          end do
        end do
      end if
      ! far: analytic 0th S' ladder kernel on the refined p nodes (target normal), fold via IPk
      do i = 1, nt
        if (cl(i)) cycle
        rt = real(tx(i),r64); zti = aimag(tx(i)); nr = real(tnx(i),r64); nz = aimag(tnx(i))
        do jp = 1, p
          rhop = real(Yp(jp),r64); zh = zti - aimag(Yp(jp)); rho = rt
          rr2 = (rho-rhop)**2 + zh*zh; chi = 1.0_r64 + rr2/(2.0_r64*rho*rhop)
          call modal_green_r64(chi, 0_8, vk, ve, Fn, An, dFn)
          Qd = -ve/(2.0_r64*(chi-1.0_r64)); SK = sqrt(rhop/rho)/twopi
          ff(jp) = SK*( -nr/(2.0_r64*rho)*vk + (nr*(rho-rhop*chi)+nz*zh)/(rho*rhop)*Qd )*wsp(jp)
        end do
        do l = 1, p
          A(i, cjo+l) = A(i, cjo+l) + sum(ff(1:p)*IPk(:,l))
        end do
      end do
    end do
    deallocate(tin, cl, ci, zc, zcn, C1, C2, C3, As, Ad, A1, A2, A3, A4, Gcq)
  end subroutine axissymlap_slpn_blockmat_r64

  subroutine axissymlap_slpn_blockmat_nmode_r64(nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, A)
    ! scalar Laplace SLPn (S'=d_n S, target-normal traction) all-modes block, real A(nt, np*p, M+1).
    ! GRADED dyadic pole refinement (the traction is ~400x more B^6-sensitive than the SLP value, so
    ! every panel obeys panlen<=0.4*pole_dist).  Close = log + smooth + TARGET Cauchy (Ad/nvp is the bare
    ! Cauchy; Sgrad=-Ad/nvp, project on the target normal zcn) p-DIRECT on the aux nodes + analytic
    ! self-node override; far = the mode-n ladder kernel.  Mirror of LapSLPnAxiBlockMat_nmode (NOT 2p).
    integer(8),   intent(in)    :: nt, p, np, M, iside, iclosed
    complex(r64), intent(in)    :: tx(nt), tnx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: sws(p*np), tpan(np+1)
    real(r64),    intent(inout) :: A(nt, np*p, M+1)
    integer(8) :: q, i, j, jp, ia, l, nc, cjo, md, k, jo, ne, npa, kk
    real(r64)  :: twopi, sumws, rho, rhop, zh, rr2, chi, rt, zti, vk, ve, Fn, An, dFn, SK, nr, nz, Qd, rn, cv
    real(r64)  :: tfirst, tlast, t0p, t1p, gam, denom, rlo, rhi, tgi, split, ff(p)
    logical    :: changed
    complex(r64) :: ic, za, zb
    real(r64)    :: tglp(p), wglp(p), Dp(p,p), IPk(p,p), IPe2(2,p), rk(p), re2(2), wsp(p)
    complex(r64) :: Yp(p), dYp(p), nvp(p), wxpp(p), ec2(2), xc(p,np)
    real(r64),    allocatable :: tin(:)
    logical,      allocatable :: cl(:)
    integer(8),   allocatable :: ci(:)
    complex(r64), allocatable :: zc(:), zcn(:), Ad(:,:)
    real(r64),    allocatable :: C1(:,:,:), C2(:,:,:), C3(:,:,:), As(:,:), A1(:,:), A2(:,:), A3(:,:), A4(:,:), Gcq(:,:)
    q = 2*p; ic = (0.0_r64,1.0_r64); twopi = 2.0_r64*acos(-1.0_r64)
    call gauss_r64(p, tglp, wglp, Dp)
    do j = 1, np
      do i = 1, p
        xc(i,j) = sx((j-1)*p+i)
      end do
    end do
    ! ---- graded dyadic pole refinement of tpan -> tin (gam=0.4 on every panel) ----
    tfirst = tpan(1)  + (1.0_r64+tglp(1))/2.0_r64*(tpan(2)-tpan(1))
    tlast  = tpan(np) + (1.0_r64+tglp(p))/2.0_r64*(tpan(np+1)-tpan(np))
    t0p = tpan(1); t1p = tpan(np+1)
    allocate(tin(np+1+1024)); ne = np+1; tin(1:ne) = tpan(1:np+1)
    do while (tin(2) > tfirst)
      split = 0.5_r64*(tin(1)+tin(2))
      do i = ne, 2, -1; tin(i+1) = tin(i); end do
      tin(2) = split; ne = ne+1
    end do
    do while (tin(ne-1) < tlast)
      split = 0.5_r64*(tin(ne)+tin(ne-1))
      tin(ne+1) = tin(ne); tin(ne) = split; ne = ne+1
    end do
    gam = 0.4_r64; changed = .true.
    do while (changed)
      changed = .false.; i = 1
      do while (i < ne)
        denom = min(tin(i)-t0p, t1p-tin(i+1))
        if (denom > 0.0_r64 .and. (tin(i+1)-tin(i)) > gam*denom) then
          do l = ne, i+1, -1; tin(l+1) = tin(l); end do
          tin(i+1) = 0.5_r64*(tin(i)+tin(i+1)); ne = ne+1; changed = .true.; i = i+2
        else
          i = i+1
        end if
      end do
    end do
    npa = ne-1
    allocate(cl(nt), ci(nt), zc(nt), zcn(nt), C1(nt,q,M+1), C2(nt,q,M+1), C3(nt,q,M+1))
    allocate(As(q,nt), Ad(q,nt), A1(q,nt), A2(q,nt), A3(q,nt), A4(q,nt), Gcq(nt,q))
    A = 0.0_r64
    do k = 1, npa
      tgi = 0.5_r64*(tin(k)+tin(k+1)); jo = 1
      do j = 1, np
        if (tgi >= tpan(j) .and. tgi < tpan(j+1)) jo = j
      end do
      denom = tpan(jo+1) - tpan(jo)
      do i = 1, p
        rlo = tin(k) + (1.0_r64+tglp(i))/2.0_r64*(tin(k+1)-tin(k))
        rk(i) = (2.0_r64*rlo - (tpan(jo)+tpan(jo+1)))/denom
      end do
      call lagrange_interp_r64(p, tglp, p, rk, IPk)
      Yp = matmul(IPk, xc(:,jo))
      rlo = (2.0_r64*tin(k)   - (tpan(jo)+tpan(jo+1)))/denom
      rhi = (2.0_r64*tin(k+1) - (tpan(jo)+tpan(jo+1)))/denom
      re2(1) = rlo; re2(2) = rhi
      call lagrange_interp_r64(p, tglp, 2_8, re2, IPe2)
      ec2 = matmul(IPe2, xc(:,jo)); za = ec2(1); zb = ec2(2)
      cjo = (jo-1)*p
      dYp = matmul(Dp, Yp)
      do jp = 1, p
        nvp(jp) = -ic*dYp(jp)/abs(dYp(jp)); wsp(jp) = wglp(jp)*abs(dYp(jp)); wxpp(jp) = dYp(jp)*wglp(jp)
      end do
      sumws = sum(wsp)
      nc = 0
      do i = 1, nt
        cl(i) = (abs(tx(i)-za) + abs(tx(i)-zb)) < 1.5_r64*sumws
        if (cl(i)) then; nc = nc + 1; ci(nc) = i; zc(nc) = tx(i); zcn(nc) = tnx(i); end if
      end do
      ! near: p-DIRECT all-mode S' close on the aux nodes (coef C1/C2/C3 + sdspecialquad As, target Cauchy
      ! -Ad/nvp) + analytic self-node override at the coincident pair, fold to coarse via IPk (== reference)
      if (nc > 0) then
        call axissymlap_slpn_coef_nmode_r64(nc, zc(1:nc), zcn(1:nc), p, Yp, M, &
             C1(1:nc,1:p,:), C2(1:nc,1:p,:), C3(1:nc,1:p,:))
        call sdspecialquad_r64(nc, zc(1:nc), p, Yp, nvp, wxpp, za, zb, iside, &
             As(1:p,1:nc), Ad(1:p,1:nc), A1(1:p,1:nc), A2(1:p,1:nc), A3(1:p,1:nc), A4(1:p,1:nc))
        do md = 0, M
          rn = 2.5_r64*log(2.0_r64)
          do kk = 1, md; rn = rn - 2.0_r64/(2.0_r64*real(kk,r64)-1.0_r64); end do
          do jp = 1, p
            do ia = 1, nc
              cv = C2(ia,jp,md+1)
              if (abs(zc(ia)-Yp(jp)) < 1.0e-13_r64) then                 ! coincident self node: analytic c2 limit
                rho = real(zc(ia),r64); nr = real(zcn(ia),r64)
                cv = nr/(2.0_r64*twopi*rho)*(1.0_r64 - rn - 0.5_r64*log(2.0_r64*rho*rho))
              end if
              Gcq(ia,jp) = twopi*C1(ia,jp,md+1)*As(jp,ia) + cv*wsp(jp) &
                         - twopi*C3(ia,jp,md+1)*real(zcn(ia)*Ad(jp,ia)/nvp(jp), r64)   ! Sgrad = -Ad/nvp
            end do
          end do
          do ia = 1, nc
            do l = 1, p
              A(ci(ia), cjo+l, md+1) = A(ci(ia), cjo+l, md+1) + sum(Gcq(ia,1:p)*IPk(:,l))
            end do
          end do
        end do
      end if
      ! far: analytic mode-n S' ladder kernel on the refined p nodes (target normal), fold via IPk
      do i = 1, nt
        if (cl(i)) cycle
        rt = real(tx(i),r64); zti = aimag(tx(i)); nr = real(tnx(i),r64); nz = aimag(tnx(i))
        block
          real(r64) :: vka(0:M), vea(0:M), Fna(0:M), Ana(0:M), dFna(0:M), vkmat(p,0:M), vemat(p,0:M)
          do jp = 1, p                                          ! carrier_all ONCE per node (was M+1 modal_green_r64 calls)
            rhop = real(Yp(jp),r64); zh = zti - aimag(Yp(jp)); rho = rt
            rr2 = (rho-rhop)**2 + zh*zh; chi = 1.0_r64 + rr2/(2.0_r64*rho*rhop)
            call modal_green_all_r64(chi, M, vka, vea, Fna, Ana, dFna)
            vkmat(jp,0:M) = vka(0:M); vemat(jp,0:M) = vea(0:M)
          end do
        do md = 0, M
          do jp = 1, p
            rhop = real(Yp(jp),r64); zh = zti - aimag(Yp(jp)); rho = rt
            rr2 = (rho-rhop)**2 + zh*zh; chi = 1.0_r64 + rr2/(2.0_r64*rho*rhop)
            vk = vkmat(jp,md); ve = vemat(jp,md)
            Qd = -ve/(2.0_r64*(chi-1.0_r64)); SK = sqrt(rhop/rho)/twopi
            ff(jp) = SK*( -nr/(2.0_r64*rho)*vk + (nr*(rho-rhop*chi)+nz*zh)/(rho*rhop)*Qd )*wsp(jp)
          end do
          do l = 1, p
            A(i, cjo+l, md+1) = A(i, cjo+l, md+1) + sum(ff(1:p)*IPk(:,l))
          end do
        end do
        end block
      end do
    end do
    deallocate(tin, cl, ci, zc, zcn, C1, C2, C3, As, Ad, A1, A2, A3, A4, Gcq)
  end subroutine axissymlap_slpn_blockmat_nmode_r64

  subroutine axissymlap_dlp_blockmat_r64(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, A)
    ! scalar Laplace DLP (double-layer potential, D=d_n' G, source-normal traction) block (mode 0), real A(nt, np*p).
    integer(8),   intent(in)    :: nt, p, np, iside, iclosed
    complex(r64), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: sws(p*np), tpan(np+1)
    real(r64),    intent(inout) :: A(nt, np*p)
    integer(8) :: q, i, j, jp, ia, l, nc, cjo, k, jo, ne, npa
    real(r64)  :: twopi, sumws, rho, rhop, zh, rr2, chi, rt, zti, vk, ve, Fn, An, dFn, SK, nr, nz, Qd, rn, cv
    real(r64)  :: tfirst, tlast, t0p, t1p, gam, denom, rlo, rhi, tgi, split, ff(p)
    logical    :: changed
    complex(r64) :: ic, za, zb
    real(r64)    :: tglp(p), wglp(p), Dp(p,p), IPk(p,p), IPe2(2,p), rk(p), re2(2), wsp(p)
    complex(r64) :: Yp(p), dYp(p), nvp(p), wxpp(p), ec2(2), xc(p,np)
    real(r64),    allocatable :: tin(:)
    logical,      allocatable :: cl(:)
    integer(8),   allocatable :: ci(:)
    complex(r64), allocatable :: zc(:), Ad(:,:)
    real(r64),    allocatable :: C1(:,:), C2(:,:), C3(:,:), As(:,:), A1(:,:), A2(:,:), A3(:,:), A4(:,:), Gcq(:,:)
    q = 2*p; ic = (0.0_r64,1.0_r64); twopi = 2.0_r64*acos(-1.0_r64)
    call gauss_r64(p, tglp, wglp, Dp)
    do j = 1, np
      do i = 1, p
        xc(i,j) = sx((j-1)*p+i)
      end do
    end do
    ! ---- graded dyadic pole refinement of tpan -> tin (gam=0.4 on every panel) ----
    tfirst = tpan(1)  + (1.0_r64+tglp(1))/2.0_r64*(tpan(2)-tpan(1))
    tlast  = tpan(np) + (1.0_r64+tglp(p))/2.0_r64*(tpan(np+1)-tpan(np))
    t0p = tpan(1); t1p = tpan(np+1)
    allocate(tin(np+1+1024)); ne = np+1; tin(1:ne) = tpan(1:np+1)
    do while (tin(2) > tfirst)
      split = 0.5_r64*(tin(1)+tin(2))
      do i = ne, 2, -1; tin(i+1) = tin(i); end do
      tin(2) = split; ne = ne+1
    end do
    do while (tin(ne-1) < tlast)
      split = 0.5_r64*(tin(ne)+tin(ne-1))
      tin(ne+1) = tin(ne); tin(ne) = split; ne = ne+1
    end do
    gam = 0.4_r64; changed = .true.
    do while (changed)
      changed = .false.; i = 1
      do while (i < ne)
        denom = min(tin(i)-t0p, t1p-tin(i+1))
        if (denom > 0.0_r64 .and. (tin(i+1)-tin(i)) > gam*denom) then
          do l = ne, i+1, -1; tin(l+1) = tin(l); end do
          tin(i+1) = 0.5_r64*(tin(i)+tin(i+1)); ne = ne+1; changed = .true.; i = i+2
        else
          i = i+1
        end if
      end do
    end do
    npa = ne-1
    allocate(cl(nt), ci(nt), zc(nt), C1(nt,q), C2(nt,q), C3(nt,q))
    allocate(As(q,nt), Ad(q,nt), A1(q,nt), A2(q,nt), A3(q,nt), A4(q,nt), Gcq(nt,q))
    A = 0.0_r64
    do k = 1, npa
      tgi = 0.5_r64*(tin(k)+tin(k+1)); jo = 1
      do j = 1, np
        if (tgi >= tpan(j) .and. tgi < tpan(j+1)) jo = j
      end do
      denom = tpan(jo+1) - tpan(jo)
      do i = 1, p
        rlo = tin(k) + (1.0_r64+tglp(i))/2.0_r64*(tin(k+1)-tin(k))
        rk(i) = (2.0_r64*rlo - (tpan(jo)+tpan(jo+1)))/denom
      end do
      call lagrange_interp_r64(p, tglp, p, rk, IPk)
      Yp = matmul(IPk, xc(:,jo))
      rlo = (2.0_r64*tin(k)   - (tpan(jo)+tpan(jo+1)))/denom
      rhi = (2.0_r64*tin(k+1) - (tpan(jo)+tpan(jo+1)))/denom
      re2(1) = rlo; re2(2) = rhi
      call lagrange_interp_r64(p, tglp, 2_8, re2, IPe2)
      ec2 = matmul(IPe2, xc(:,jo)); za = ec2(1); zb = ec2(2)
      cjo = (jo-1)*p
      dYp = matmul(Dp, Yp)
      do jp = 1, p
        nvp(jp) = -ic*dYp(jp)/abs(dYp(jp)); wsp(jp) = wglp(jp)*abs(dYp(jp)); wxpp(jp) = dYp(jp)*wglp(jp)
      end do
      sumws = sum(wsp)
      nc = 0
      do i = 1, nt
        cl(i) = (abs(tx(i)-za) + abs(tx(i)-zb)) < 1.5_r64*sumws
        if (cl(i)) then; nc = nc + 1; ci(nc) = i; zc(nc) = tx(i); end if
      end do
      ! near: p-DIRECT 0th DLP close (coef + sdspecialquad As, SOURCE Cauchy +2pi c1S real(Ad)) + override, fold via IPk
      if (nc > 0) then
        call axissymlap_dlp_coef_r64(nc, zc(1:nc), p, Yp, nvp, C1(1:nc,1:p), C2(1:nc,1:p), C3(1:nc,1:p))
        call sdspecialquad_r64(nc, zc(1:nc), p, Yp, nvp, wxpp, za, zb, iside, &
             As(1:p,1:nc), Ad(1:p,1:nc), A1(1:p,1:nc), A2(1:p,1:nc), A3(1:p,1:nc), A4(1:p,1:nc))
        rn = 2.5_r64*log(2.0_r64)
        do jp = 1, p
          do ia = 1, nc
            cv = C2(ia,jp)
            if (abs(zc(ia)-Yp(jp)) < 1.0e-13_r64) then                   ! coincident self node: analytic c2 limit
              rhop = real(Yp(jp),r64); nr = real(nvp(jp),r64)
              cv = nr/(2.0_r64*twopi*rhop)*(1.0_r64 - rn - 0.5_r64*log(2.0_r64*rhop*rhop))
            end if
            Gcq(ia,jp) = twopi*C1(ia,jp)*As(jp,ia) + cv*wsp(jp) &
                       + twopi*C3(ia,jp)*real(Ad(jp,ia), r64)             ! source Cauchy = +2pi c1S real(Dval)
          end do
        end do
        do ia = 1, nc
          do l = 1, p
            A(ci(ia), cjo+l) = A(ci(ia), cjo+l) + sum(Gcq(ia,1:p)*IPk(:,l))
          end do
        end do
      end if
      ! far: analytic 0th DLP ladder kernel on the aux nodes (source normal nvp), fold via IPk
      do i = 1, nt
        if (cl(i)) cycle
        rt = real(tx(i),r64); zti = aimag(tx(i))
        do jp = 1, p
          rhop = real(Yp(jp),r64); zh = zti - aimag(Yp(jp)); rho = rt
          nr = real(nvp(jp),r64); nz = aimag(nvp(jp))
          rr2 = (rho-rhop)**2 + zh*zh; chi = 1.0_r64 + rr2/(2.0_r64*rho*rhop)
          call modal_green_r64(chi, 0_8, vk, ve, Fn, An, dFn)
          Qd = -ve/(2.0_r64*(chi-1.0_r64)); SK = sqrt(rhop/rho)/twopi
          ff(jp) = SK*( -nr/(2.0_r64*rhop)*vk + (nr*(rhop-rho*chi)-nz*zh)/(rho*rhop)*Qd )*wsp(jp)
        end do
        do l = 1, p
          A(i, cjo+l) = A(i, cjo+l) + sum(ff(1:p)*IPk(:,l))
        end do
      end do
    end do
    deallocate(tin, cl, ci, zc, C1, C2, C3, As, Ad, A1, A2, A3, A4, Gcq)
  end subroutine axissymlap_dlp_blockmat_r64

  subroutine axissymlap_dlp_blockmat_nmode_r64(nt, tx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, A)
    ! scalar Laplace DLP (double-layer potential, D=d_n' G, source-normal traction) all-modes block,
    ! real A(nt, np*p, M+1).  Source-normal mirror of the SLPn block: GRADED dyadic pole refinement +
    ! p-DIRECT close on the aux nodes (log + smooth + SOURCE Cauchy +2pi c1S real(Ad)) + analytic self-node
    ! override; far = the mode-n ladder kernel with the source normal.  Mirror of LapDLPAxiBlockMat_nmode.
    integer(8),   intent(in)    :: nt, p, np, M, iside, iclosed
    complex(r64), intent(in)    :: tx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: sws(p*np), tpan(np+1)
    real(r64),    intent(inout) :: A(nt, np*p, M+1)
    integer(8) :: q, i, j, jp, ia, l, nc, cjo, md, k, jo, ne, npa, kk
    real(r64)  :: twopi, sumws, rho, rhop, zh, rr2, chi, rt, zti, vk, ve, Fn, An, dFn, SK, nr, nz, Qd, rn, cv
    real(r64)  :: tfirst, tlast, t0p, t1p, gam, denom, rlo, rhi, tgi, split, ff(p)
    logical    :: changed
    complex(r64) :: ic, za, zb
    real(r64)    :: tglp(p), wglp(p), Dp(p,p), IPk(p,p), IPe2(2,p), rk(p), re2(2), wsp(p)
    complex(r64) :: Yp(p), dYp(p), nvp(p), wxpp(p), ec2(2), xc(p,np)
    real(r64),    allocatable :: tin(:)
    logical,      allocatable :: cl(:)
    integer(8),   allocatable :: ci(:)
    complex(r64), allocatable :: zc(:), Ad(:,:)
    real(r64),    allocatable :: C1(:,:,:), C2(:,:,:), C3(:,:,:), As(:,:), A1(:,:), A2(:,:), A3(:,:), A4(:,:), Gcq(:,:)
    q = 2*p; ic = (0.0_r64,1.0_r64); twopi = 2.0_r64*acos(-1.0_r64)
    call gauss_r64(p, tglp, wglp, Dp)
    do j = 1, np
      do i = 1, p
        xc(i,j) = sx((j-1)*p+i)
      end do
    end do
    ! ---- graded dyadic pole refinement of tpan -> tin (gam=0.4 on every panel) ----
    tfirst = tpan(1)  + (1.0_r64+tglp(1))/2.0_r64*(tpan(2)-tpan(1))
    tlast  = tpan(np) + (1.0_r64+tglp(p))/2.0_r64*(tpan(np+1)-tpan(np))
    t0p = tpan(1); t1p = tpan(np+1)
    allocate(tin(np+1+1024)); ne = np+1; tin(1:ne) = tpan(1:np+1)
    do while (tin(2) > tfirst)
      split = 0.5_r64*(tin(1)+tin(2))
      do i = ne, 2, -1; tin(i+1) = tin(i); end do
      tin(2) = split; ne = ne+1
    end do
    do while (tin(ne-1) < tlast)
      split = 0.5_r64*(tin(ne)+tin(ne-1))
      tin(ne+1) = tin(ne); tin(ne) = split; ne = ne+1
    end do
    gam = 0.4_r64; changed = .true.
    do while (changed)
      changed = .false.; i = 1
      do while (i < ne)
        denom = min(tin(i)-t0p, t1p-tin(i+1))
        if (denom > 0.0_r64 .and. (tin(i+1)-tin(i)) > gam*denom) then
          do l = ne, i+1, -1; tin(l+1) = tin(l); end do
          tin(i+1) = 0.5_r64*(tin(i)+tin(i+1)); ne = ne+1; changed = .true.; i = i+2
        else
          i = i+1
        end if
      end do
    end do
    npa = ne-1
    allocate(cl(nt), ci(nt), zc(nt), C1(nt,q,M+1), C2(nt,q,M+1), C3(nt,q,M+1))
    allocate(As(q,nt), Ad(q,nt), A1(q,nt), A2(q,nt), A3(q,nt), A4(q,nt), Gcq(nt,q))
    A = 0.0_r64
    do k = 1, npa
      tgi = 0.5_r64*(tin(k)+tin(k+1)); jo = 1
      do j = 1, np
        if (tgi >= tpan(j) .and. tgi < tpan(j+1)) jo = j
      end do
      denom = tpan(jo+1) - tpan(jo)
      do i = 1, p
        rlo = tin(k) + (1.0_r64+tglp(i))/2.0_r64*(tin(k+1)-tin(k))
        rk(i) = (2.0_r64*rlo - (tpan(jo)+tpan(jo+1)))/denom
      end do
      call lagrange_interp_r64(p, tglp, p, rk, IPk)
      Yp = matmul(IPk, xc(:,jo))
      rlo = (2.0_r64*tin(k)   - (tpan(jo)+tpan(jo+1)))/denom
      rhi = (2.0_r64*tin(k+1) - (tpan(jo)+tpan(jo+1)))/denom
      re2(1) = rlo; re2(2) = rhi
      call lagrange_interp_r64(p, tglp, 2_8, re2, IPe2)
      ec2 = matmul(IPe2, xc(:,jo)); za = ec2(1); zb = ec2(2)
      cjo = (jo-1)*p
      dYp = matmul(Dp, Yp)
      do jp = 1, p
        nvp(jp) = -ic*dYp(jp)/abs(dYp(jp)); wsp(jp) = wglp(jp)*abs(dYp(jp)); wxpp(jp) = dYp(jp)*wglp(jp)
      end do
      sumws = sum(wsp)
      nc = 0
      do i = 1, nt
        cl(i) = (abs(tx(i)-za) + abs(tx(i)-zb)) < 1.5_r64*sumws
        if (cl(i)) then; nc = nc + 1; ci(nc) = i; zc(nc) = tx(i); end if
      end do
      ! near: p-DIRECT all-mode DLP close on the aux nodes (coef C1/C2/C3 + sdspecialquad As, SOURCE Cauchy
      ! +2pi c1S real(Ad)) + analytic self-node override, fold to coarse via IPk (== LapDLPAxiBlockMat_nmode)
      if (nc > 0) then
        call axissymlap_dlp_coef_nmode_r64(nc, zc(1:nc), p, Yp, nvp, M, &
             C1(1:nc,1:p,:), C2(1:nc,1:p,:), C3(1:nc,1:p,:))
        call sdspecialquad_r64(nc, zc(1:nc), p, Yp, nvp, wxpp, za, zb, iside, &
             As(1:p,1:nc), Ad(1:p,1:nc), A1(1:p,1:nc), A2(1:p,1:nc), A3(1:p,1:nc), A4(1:p,1:nc))
        do md = 0, M
          rn = 2.5_r64*log(2.0_r64)
          do kk = 1, md; rn = rn - 2.0_r64/(2.0_r64*real(kk,r64)-1.0_r64); end do
          do jp = 1, p
            do ia = 1, nc
              cv = C2(ia,jp,md+1)
              if (abs(zc(ia)-Yp(jp)) < 1.0e-13_r64) then                 ! coincident self node: analytic c2 limit
                rhop = real(Yp(jp),r64); nr = real(nvp(jp),r64)
                cv = nr/(2.0_r64*twopi*rhop)*(1.0_r64 - rn - 0.5_r64*log(2.0_r64*rhop*rhop))
              end if
              Gcq(ia,jp) = twopi*C1(ia,jp,md+1)*As(jp,ia) + cv*wsp(jp) &
                         + twopi*C3(ia,jp,md+1)*real(Ad(jp,ia), r64)      ! source Cauchy = +2pi c1S real(Dval)
            end do
          end do
          do ia = 1, nc
            do l = 1, p
              A(ci(ia), cjo+l, md+1) = A(ci(ia), cjo+l, md+1) + sum(Gcq(ia,1:p)*IPk(:,l))
            end do
          end do
        end do
      end if
      ! far: analytic mode-n DLP ladder kernel on the aux nodes (source normal nvp), fold via IPk
      do i = 1, nt
        if (cl(i)) cycle
        rt = real(tx(i),r64); zti = aimag(tx(i))
        block
          real(r64) :: vka(0:M), vea(0:M), Fna(0:M), Ana(0:M), dFna(0:M), vkmat(p,0:M), vemat(p,0:M)
          do jp = 1, p                                          ! carrier_all ONCE per node (was M+1 modal_green_r64 calls)
            rhop = real(Yp(jp),r64); zh = zti - aimag(Yp(jp)); rho = rt
            rr2 = (rho-rhop)**2 + zh*zh; chi = 1.0_r64 + rr2/(2.0_r64*rho*rhop)
            call modal_green_all_r64(chi, M, vka, vea, Fna, Ana, dFna)
            vkmat(jp,0:M) = vka(0:M); vemat(jp,0:M) = vea(0:M)
          end do
        do md = 0, M
          do jp = 1, p
            rhop = real(Yp(jp),r64); zh = zti - aimag(Yp(jp)); rho = rt
            nr = real(nvp(jp),r64); nz = aimag(nvp(jp))
            rr2 = (rho-rhop)**2 + zh*zh; chi = 1.0_r64 + rr2/(2.0_r64*rho*rhop)
            vk = vkmat(jp,md); ve = vemat(jp,md)
            Qd = -ve/(2.0_r64*(chi-1.0_r64)); SK = sqrt(rhop/rho)/twopi
            ff(jp) = SK*( -nr/(2.0_r64*rhop)*vk + (nr*(rhop-rho*chi)-nz*zh)/(rho*rhop)*Qd )*wsp(jp)
          end do
          do l = 1, p
            A(i, cjo+l, md+1) = A(i, cjo+l, md+1) + sum(ff(1:p)*IPk(:,l))
          end do
        end do
        end block
      end do
    end do
    deallocate(tin, cl, ci, zc, C1, C2, C3, As, Ad, A1, A2, A3, A4, Gcq)
  end subroutine axissymlap_dlp_blockmat_nmode_r64

  subroutine axissymlap_dlpn_blockmat_r64(nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, iside, iclosed, A)
    ! scalar Laplace DLPn (D'=d_n d_n' G, hypersingular double-layer traction) block (mode 0), real A(nt, np*p).
    integer(8),   intent(in)    :: nt, p, np, iside, iclosed
    complex(r64), intent(in)    :: tx(nt), tnx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: sws(p*np), tpan(np+1)
    real(r64),    intent(inout) :: A(nt, np*p)
    integer(8) :: q, i, j, jp, iq, ia, l, nc, cjo, k, jo, ne, npa
    real(r64)  :: twopi, sumws, rho, rhop, zh, drho, rr2, tchi, chi, rt, zti, vk, ve, Fn, An, dFn, Qd, Qdd
    real(r64)  :: nr, nz, nrp, nzp, Bq, Br, Brp, Brrp, chir, chiz, chirp, chizp, chirrp, chirzp, chizrp, chizzp
    real(r64)  :: Grrp, Grzp, Gzrp, Gzzp, rowsum
    real(r64)  :: tfirst, tlast, t0p, t1p, gam, denom, rlo, rhi, tgi, split, ff(p)
    logical    :: changed, isself
    complex(r64) :: ic, za, zb, Azc
    real(r64)    :: tglp(p), wglp(p), Dp(p,p), tglq(2*p), wglq(2*p), Dq(2*p,2*p), IP2(2*p,p)
    real(r64)    :: IPk(p,p), IPe2(2,p), rk(p), re2(2), wsq(2*p), wsp(p)
    complex(r64) :: Yp(p), Yq(2*p), dYp(p), dYq(2*p), nvp(p), nvq(2*p), wxpq(2*p), ec2(2), xc(p,np)
    real(r64),    allocatable :: tin(:)
    logical,      allocatable :: cl(:)
    integer(8),   allocatable :: ci(:)
    complex(r64), allocatable :: zc(:), zcn(:), Ad(:,:)
    real(r64),    allocatable :: C1(:,:), C2(:,:), C3a(:,:), C3b(:,:), C4(:,:)
    real(r64),    allocatable :: As(:,:), A1(:,:), A2(:,:), A3(:,:), A4(:,:), Gcq(:,:), Gpc(:,:)
    q = 2*p; ic = (0.0_r64,1.0_r64); twopi = 2.0_r64*acos(-1.0_r64)
    call gauss_r64(p, tglp, wglp, Dp)
    call gauss_r64(2_8*p, tglq, wglq, Dq)
    call lagrange_interp_r64(p, tglp, q, tglq, IP2)
    do j = 1, np
      do i = 1, p
        xc(i,j) = sx((j-1)*p+i)
      end do
    end do
    isself = (nt == np*p)
    if (isself) then
      do i = 1, nt
        if (abs(tx(i)-sx(i)) > 1.0e-13_r64) isself = .false.
      end do
    end if
    tfirst = tpan(1)  + (1.0_r64+tglp(1))/2.0_r64*(tpan(2)-tpan(1))
    tlast  = tpan(np) + (1.0_r64+tglp(p))/2.0_r64*(tpan(np+1)-tpan(np))
    t0p = tpan(1); t1p = tpan(np+1)
    allocate(tin(np+1+1024)); ne = np+1; tin(1:ne) = tpan(1:np+1)
    do while (tin(2) > tfirst)
      split = 0.5_r64*(tin(1)+tin(2))
      do i = ne, 2, -1; tin(i+1) = tin(i); end do
      tin(2) = split; ne = ne+1
    end do
    do while (tin(ne-1) < tlast)
      split = 0.5_r64*(tin(ne)+tin(ne-1))
      tin(ne+1) = tin(ne); tin(ne) = split; ne = ne+1
    end do
    gam = 0.4_r64; changed = .true.
    do while (changed)
      changed = .false.; i = 1
      do while (i < ne)
        denom = min(tin(i)-t0p, t1p-tin(i+1))
        if (denom > 0.0_r64 .and. (tin(i+1)-tin(i)) > gam*denom) then
          do l = ne, i+1, -1; tin(l+1) = tin(l); end do
          tin(i+1) = 0.5_r64*(tin(i)+tin(i+1)); ne = ne+1; changed = .true.; i = i+2
        else
          i = i+1
        end if
      end do
    end do
    npa = ne-1
    allocate(cl(nt), ci(nt), zc(nt), zcn(nt), C1(nt,q), C2(nt,q), C3a(nt,q), C3b(nt,q), C4(nt,q))
    allocate(As(q,nt), Ad(q,nt), A1(q,nt), A2(q,nt), A3(q,nt), A4(q,nt), Gcq(nt,q), Gpc(nt,p))
    A = 0.0_r64
    do k = 1, npa
      tgi = 0.5_r64*(tin(k)+tin(k+1)); jo = 1
      do j = 1, np
        if (tgi >= tpan(j) .and. tgi < tpan(j+1)) jo = j
      end do
      denom = tpan(jo+1) - tpan(jo)
      do i = 1, p
        rlo = tin(k) + (1.0_r64+tglp(i))/2.0_r64*(tin(k+1)-tin(k))
        rk(i) = (2.0_r64*rlo - (tpan(jo)+tpan(jo+1)))/denom
      end do
      call lagrange_interp_r64(p, tglp, p, rk, IPk)
      Yp = matmul(IPk, xc(:,jo))
      rlo = (2.0_r64*tin(k)   - (tpan(jo)+tpan(jo+1)))/denom
      rhi = (2.0_r64*tin(k+1) - (tpan(jo)+tpan(jo+1)))/denom
      re2(1) = rlo; re2(2) = rhi
      call lagrange_interp_r64(p, tglp, 2_8, re2, IPe2)
      ec2 = matmul(IPe2, xc(:,jo)); za = ec2(1); zb = ec2(2)
      cjo = (jo-1)*p
      Yq = matmul(IP2, Yp); dYq = matmul(Dq, Yq)
      do iq = 1, q
        nvq(iq) = -ic*dYq(iq)/abs(dYq(iq)); wsq(iq) = wglq(iq)*abs(dYq(iq)); wxpq(iq) = dYq(iq)*wglq(iq)
      end do
      dYp = matmul(Dp, Yp)
      do jp = 1, p
        nvp(jp) = -ic*dYp(jp)/abs(dYp(jp)); wsp(jp) = wglp(jp)*abs(dYp(jp))
      end do
      sumws = sum(wsp)
      nc = 0
      do i = 1, nt
        cl(i) = (abs(tx(i)-za) + abs(tx(i)-zb)) < 1.5_r64*sumws
        if (cl(i)) then; nc = nc + 1; ci(nc) = i; zc(nc) = tx(i); zcn(nc) = tnx(i); end if
      end do
      ! near: 2p UPSAMPLE 0th DLPn close (coef 5-bucket + sdspecialquad: log As, source Cauchy Ad, deriv Az),
      ! project 2p->p (IP2), fold via IPk  (== LapDLPnAxiBlockMat_nmode)
      if (nc > 0) then
        call axissymlap_dlpn_coef_r64(nc, zc(1:nc), zcn(1:nc), q, Yq, nvq, &
             C1(1:nc,1:q), C2(1:nc,1:q), C3a(1:nc,1:q), C3b(1:nc,1:q), C4(1:nc,1:q))
        call sdspecialquad_r64(nc, zc(1:nc), q, Yq, nvq, wxpq, za, zb, iside, &
             As(1:q,1:nc), Ad(1:q,1:nc), A1(1:q,1:nc), A2(1:q,1:nc), A3(1:q,1:nc), A4(1:q,1:nc))
        do iq = 1, q
          do ia = 1, nc
            Azc = cmplx(A1(iq,ia), -A2(iq,ia), r64)
            Gcq(ia,iq) = twopi*C1(ia,iq)*As(iq,ia) + C2(ia,iq)*wsq(iq) &
                       - twopi*C3a(ia,iq)*real(zcn(ia)*Ad(iq,ia)/nvq(iq), r64) &
                       + twopi*C3b(ia,iq)*real(Ad(iq,ia), r64) &
                       + twopi*C4(ia,iq)*real(zcn(ia)*Azc, r64)
          end do
        end do
        Gpc(1:nc,1:p) = matmul(Gcq(1:nc,1:q), IP2)
        do ia = 1, nc
          do l = 1, p
            A(ci(ia), cjo+l) = A(ci(ia), cjo+l) + sum(Gpc(ia,1:p)*IPk(:,l))
          end do
        end do
      end if
      ! far: (Q,Q',Q'') Hessian contraction of G^phys=V_K/sqrt(rho rho') on the p nodes (both normals), fold via IPk
      do i = 1, nt
        if (cl(i)) cycle
        rt = real(tx(i),r64); zti = aimag(tx(i)); nr = real(tnx(i),r64); nz = aimag(tnx(i))
        do jp = 1, p
          rhop = real(Yp(jp),r64); zh = zti - aimag(Yp(jp)); rho = rt; drho = rho - rhop
          nrp = real(nvp(jp),r64); nzp = aimag(nvp(jp))
          rr2 = drho*drho + zh*zh; tchi = rr2/(2.0_r64*rho*rhop); chi = 1.0_r64 + tchi
          call modal_green_r64(chi, 0_8, vk, ve, Fn, An, dFn)
          Qd = -ve/(2.0_r64*tchi); Qdd = (2.0_r64*chi*Qd + 0.25_r64*vk)/(1.0_r64-chi*chi)
          Bq = 1.0_r64/sqrt(rho*rhop); Br = -Bq/(2.0_r64*rho); Brp = -Bq/(2.0_r64*rhop); Brrp = Bq/(4.0_r64*rho*rhop)
          chir = drho/(rho*rhop) - rr2/(2.0_r64*rho*rho*rhop)
          chiz = zh/(rho*rhop)
          chirp = -drho/(rho*rhop) - rr2/(2.0_r64*rho*rhop*rhop)
          chizp = -zh/(rho*rhop)
          chirrp = -1.0_r64/(rhop*rhop) + drho/(rho*rho*rhop) + rr2/(2.0_r64*rho*rho*rhop*rhop)
          chirzp = zh/(rho*rho*rhop)
          chizrp = -zh/(rho*rhop*rhop)
          chizzp = -1.0_r64/(rho*rhop)
          Grrp = Brrp*vk + Br*Qd*chirp + Brp*Qd*chir + Bq*Qdd*chir*chirp + Bq*Qd*chirrp
          Grzp = Br*Qd*chizp + Bq*Qdd*chir*chizp + Bq*Qd*chirzp
          Gzrp = Brp*Qd*chiz + Bq*Qdd*chiz*chirp + Bq*Qd*chizrp
          Gzzp = Bq*Qdd*chiz*chizp + Bq*Qd*chizzp
          ff(jp) = (rhop/twopi)*( nr*nrp*Grrp + nr*nzp*Grzp + nz*nrp*Gzrp + nz*nzp*Gzzp )*wsp(jp)
        end do
        do l = 1, p
          A(i, cjo+l) = A(i, cjo+l) + sum(ff(1:p)*IPk(:,l))
        end do
      end do
    end do
    if (isself) then                                              ! m=0 self: D'_0[1] = 0 (hypersingular nullspace)
      do i = 1, nt
        rowsum = sum(A(i,1:np*p))
        A(i,i) = A(i,i) - rowsum
      end do
    end if
    deallocate(tin, cl, ci, zc, zcn, C1, C2, C3a, C3b, C4, As, Ad, A1, A2, A3, A4, Gcq, Gpc)
  end subroutine axissymlap_dlpn_blockmat_r64

  subroutine axissymlap_dlpn_blockmat_nmode_r64(nt, tx, tnx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, A)
    ! scalar Laplace DLPn (D'=d_n d_n' G, hypersingular double-layer traction) all-modes block, real
    ! A(nt, np*p, M+1).  Both normals (target tnx, source snx).  GRADED dyadic pole refinement + 2p
    ! UPSAMPLE-AND-PROJECT close (a p target never coincides with a 2p source) with the 5-bucket split
    ! (log + smooth + target Cauchy -Ad/nvq + source Cauchy +Ad + deriv-Cauchy Az); far = the (Q,Q',Q'')
    ! Hessian contraction of G^phys=V_K/sqrt(rho rho').  m=0 SELF row-sum fix D'_0[1]=0.  Mirror of
    ! LapDLPnAxiBlockMat_nmode.
    integer(8),   intent(in)    :: nt, p, np, M, iside, iclosed
    complex(r64), intent(in)    :: tx(nt), tnx(nt), sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: sws(p*np), tpan(np+1)
    real(r64),    intent(inout) :: A(nt, np*p, M+1)
    integer(8) :: q, i, j, jp, iq, ia, l, nc, cjo, md, k, jo, ne, npa
    real(r64)  :: twopi, sumws, rho, rhop, zh, drho, rr2, tchi, chi, rt, zti, vk, ve, Fn, An, dFn, Qd, Qdd
    real(r64)  :: nr, nz, nrp, nzp, Bq, Br, Brp, Brrp, chir, chiz, chirp, chizp, chirrp, chirzp, chizrp, chizzp
    real(r64)  :: Grrp, Grzp, Gzrp, Gzzp, rowsum, rn2
    real(r64)  :: tfirst, tlast, t0p, t1p, gam, denom, rlo, rhi, tgi, split, ff(p)
    logical    :: changed, isself
    complex(r64) :: ic, za, zb, Azc
    real(r64)    :: tglp(p), wglp(p), Dp(p,p), tglq(2*p), wglq(2*p), Dq(2*p,2*p), IP2(2*p,p)
    real(r64)    :: IPk(p,p), IPe2(2,p), rk(p), re2(2), wsq(2*p), wsp(p)
    complex(r64) :: Yp(p), Yq(2*p), dYp(p), dYq(2*p), nvp(p), nvq(2*p), wxpq(2*p), ec2(2), xc(p,np)
    real(r64),    allocatable :: tin(:)
    logical,      allocatable :: cl(:)
    integer(8),   allocatable :: ci(:)
    complex(r64), allocatable :: zc(:), zcn(:), Ad(:,:)
    real(r64),    allocatable :: C1(:,:,:), C2(:,:,:), C3a(:,:,:), C3b(:,:,:), C4(:,:,:)
    real(r64),    allocatable :: As(:,:), A1(:,:), A2(:,:), A3(:,:), A4(:,:), Gcq(:,:), Gpc(:,:)
    q = 2*p; ic = (0.0_r64,1.0_r64); twopi = 2.0_r64*acos(-1.0_r64)
    call gauss_r64(p, tglp, wglp, Dp)
    call gauss_r64(2_8*p, tglq, wglq, Dq)
    call lagrange_interp_r64(p, tglp, q, tglq, IP2)
    do j = 1, np
      do i = 1, p
        xc(i,j) = sx((j-1)*p+i)
      end do
    end do
    isself = (nt == np*p)
    if (isself) then
      do i = 1, nt
        if (abs(tx(i)-sx(i)) > 1.0e-13_r64) isself = .false.
      end do
    end if
    tfirst = tpan(1)  + (1.0_r64+tglp(1))/2.0_r64*(tpan(2)-tpan(1))
    tlast  = tpan(np) + (1.0_r64+tglp(p))/2.0_r64*(tpan(np+1)-tpan(np))
    t0p = tpan(1); t1p = tpan(np+1)
    allocate(tin(np+1+1024)); ne = np+1; tin(1:ne) = tpan(1:np+1)
    do while (tin(2) > tfirst)
      split = 0.5_r64*(tin(1)+tin(2))
      do i = ne, 2, -1; tin(i+1) = tin(i); end do
      tin(2) = split; ne = ne+1
    end do
    do while (tin(ne-1) < tlast)
      split = 0.5_r64*(tin(ne)+tin(ne-1))
      tin(ne+1) = tin(ne); tin(ne) = split; ne = ne+1
    end do
    gam = 0.4_r64; changed = .true.
    do while (changed)
      changed = .false.; i = 1
      do while (i < ne)
        denom = min(tin(i)-t0p, t1p-tin(i+1))
        if (denom > 0.0_r64 .and. (tin(i+1)-tin(i)) > gam*denom) then
          do l = ne, i+1, -1; tin(l+1) = tin(l); end do
          tin(i+1) = 0.5_r64*(tin(i)+tin(i+1)); ne = ne+1; changed = .true.; i = i+2
        else
          i = i+1
        end if
      end do
    end do
    npa = ne-1
    allocate(cl(nt), ci(nt), zc(nt), zcn(nt))
    allocate(C1(nt,q,M+1), C2(nt,q,M+1), C3a(nt,q,M+1), C3b(nt,q,M+1), C4(nt,q,M+1))
    allocate(As(q,nt), Ad(q,nt), A1(q,nt), A2(q,nt), A3(q,nt), A4(q,nt), Gcq(nt,q), Gpc(nt,p))
    A = 0.0_r64
    do k = 1, npa
      tgi = 0.5_r64*(tin(k)+tin(k+1)); jo = 1
      do j = 1, np
        if (tgi >= tpan(j) .and. tgi < tpan(j+1)) jo = j
      end do
      denom = tpan(jo+1) - tpan(jo)
      do i = 1, p
        rlo = tin(k) + (1.0_r64+tglp(i))/2.0_r64*(tin(k+1)-tin(k))
        rk(i) = (2.0_r64*rlo - (tpan(jo)+tpan(jo+1)))/denom
      end do
      call lagrange_interp_r64(p, tglp, p, rk, IPk)
      Yp = matmul(IPk, xc(:,jo))
      rlo = (2.0_r64*tin(k)   - (tpan(jo)+tpan(jo+1)))/denom
      rhi = (2.0_r64*tin(k+1) - (tpan(jo)+tpan(jo+1)))/denom
      re2(1) = rlo; re2(2) = rhi
      call lagrange_interp_r64(p, tglp, 2_8, re2, IPe2)
      ec2 = matmul(IPe2, xc(:,jo)); za = ec2(1); zb = ec2(2)
      cjo = (jo-1)*p
      Yq = matmul(IP2, Yp); dYq = matmul(Dq, Yq)
      do iq = 1, q
        nvq(iq) = -ic*dYq(iq)/abs(dYq(iq)); wsq(iq) = wglq(iq)*abs(dYq(iq)); wxpq(iq) = dYq(iq)*wglq(iq)
      end do
      dYp = matmul(Dp, Yp)
      do jp = 1, p
        nvp(jp) = -ic*dYp(jp)/abs(dYp(jp)); wsp(jp) = wglp(jp)*abs(dYp(jp))
      end do
      sumws = sum(wsp)
      nc = 0
      do i = 1, nt
        cl(i) = (abs(tx(i)-za) + abs(tx(i)-zb)) < 1.5_r64*sumws
        if (cl(i)) then; nc = nc + 1; ci(nc) = i; zc(nc) = tx(i); zcn(nc) = tnx(i); end if
      end do
      ! near: 2p UPSAMPLE all-mode DLPn close (coef 5-bucket + sdspecialquad), project 2p->p (IP2), fold via IPk
      if (nc > 0) then
        call axissymlap_dlpn_coef_nmode_r64(nc, zc(1:nc), zcn(1:nc), q, Yq, nvq, M, &
             C1(1:nc,1:q,:), C2(1:nc,1:q,:), C3a(1:nc,1:q,:), C3b(1:nc,1:q,:), C4(1:nc,1:q,:))
        call sdspecialquad_r64(nc, zc(1:nc), q, Yq, nvq, wxpq, za, zb, iside, &
             As(1:q,1:nc), Ad(1:q,1:nc), A1(1:q,1:nc), A2(1:q,1:nc), A3(1:q,1:nc), A4(1:q,1:nc))
        do md = 0, M
          do iq = 1, q
            do ia = 1, nc
              Azc = cmplx(A1(iq,ia), -A2(iq,ia), r64)
              Gcq(ia,iq) = twopi*C1(ia,iq,md+1)*As(iq,ia) + C2(ia,iq,md+1)*wsq(iq) &
                         - twopi*C3a(ia,iq,md+1)*real(zcn(ia)*Ad(iq,ia)/nvq(iq), r64) &
                         + twopi*C3b(ia,iq,md+1)*real(Ad(iq,ia), r64) &
                         + twopi*C4(ia,iq,md+1)*real(zcn(ia)*Azc, r64)
            end do
          end do
          Gpc(1:nc,1:p) = matmul(Gcq(1:nc,1:q), IP2)
          do ia = 1, nc
            do l = 1, p
              A(ci(ia), cjo+l, md+1) = A(ci(ia), cjo+l, md+1) + sum(Gpc(ia,1:p)*IPk(:,l))
            end do
          end do
        end do
      end if
      ! far: (Q,Q',Q'') Hessian contraction on the p nodes (both normals), fold via IPk
      do i = 1, nt
        if (cl(i)) cycle
        rt = real(tx(i),r64); zti = aimag(tx(i)); nr = real(tnx(i),r64); nz = aimag(tnx(i))
        block
          real(r64) :: vka(0:M), vea(0:M), Fna(0:M), Ana(0:M), dFna(0:M), vkmat(p,0:M), vemat(p,0:M)
          do jp = 1, p                                          ! carrier_all ONCE per node (was M+1 modal_green_r64 calls)
            rhop = real(Yp(jp),r64); zh = zti - aimag(Yp(jp)); rho = rt
            rr2 = (rho-rhop)**2 + zh*zh; chi = 1.0_r64 + rr2/(2.0_r64*rho*rhop)
            call modal_green_all_r64(chi, M, vka, vea, Fna, Ana, dFna)
            vkmat(jp,0:M) = vka(0:M); vemat(jp,0:M) = vea(0:M)
          end do
        do md = 0, M
          rn2 = real(md,r64)**2
          do jp = 1, p
            rhop = real(Yp(jp),r64); zh = zti - aimag(Yp(jp)); rho = rt; drho = rho - rhop
            nrp = real(nvp(jp),r64); nzp = aimag(nvp(jp))
            rr2 = drho*drho + zh*zh; tchi = rr2/(2.0_r64*rho*rhop); chi = 1.0_r64 + tchi
            vk = vkmat(jp,md); ve = vemat(jp,md)
            Qd = -ve/(2.0_r64*tchi); Qdd = (2.0_r64*chi*Qd - (rn2-0.25_r64)*vk)/(1.0_r64-chi*chi)
            Bq = 1.0_r64/sqrt(rho*rhop); Br = -Bq/(2.0_r64*rho); Brp = -Bq/(2.0_r64*rhop); Brrp = Bq/(4.0_r64*rho*rhop)
            chir = drho/(rho*rhop) - rr2/(2.0_r64*rho*rho*rhop)
            chiz = zh/(rho*rhop)
            chirp = -drho/(rho*rhop) - rr2/(2.0_r64*rho*rhop*rhop)
            chizp = -zh/(rho*rhop)
            chirrp = -1.0_r64/(rhop*rhop) + drho/(rho*rho*rhop) + rr2/(2.0_r64*rho*rho*rhop*rhop)
            chirzp = zh/(rho*rho*rhop)
            chizrp = -zh/(rho*rhop*rhop)
            chizzp = -1.0_r64/(rho*rhop)
            Grrp = Brrp*vk + Br*Qd*chirp + Brp*Qd*chir + Bq*Qdd*chir*chirp + Bq*Qd*chirrp
            Grzp = Br*Qd*chizp + Bq*Qdd*chir*chizp + Bq*Qd*chirzp
            Gzrp = Brp*Qd*chiz + Bq*Qdd*chiz*chirp + Bq*Qd*chizrp
            Gzzp = Bq*Qdd*chiz*chizp + Bq*Qd*chizzp
            ff(jp) = (rhop/twopi)*( nr*nrp*Grrp + nr*nzp*Grzp + nz*nrp*Gzrp + nz*nzp*Gzzp )*wsp(jp)
          end do
          do l = 1, p
            A(i, cjo+l, md+1) = A(i, cjo+l, md+1) + sum(ff(1:p)*IPk(:,l))
          end do
        end do
        end block
      end do
    end do
    if (isself) then                                              ! m=0 self: D'_0[1] = 0 (hypersingular nullspace)
      do i = 1, nt
        rowsum = sum(A(i,1:np*p,1))
        A(i,i,1) = A(i,i,1) - rowsum
      end do
    end if
    deallocate(tin, cl, ci, zc, zcn, C1, C2, C3a, C3b, C4, As, Ad, A1, A2, A3, A4, Gcq, Gpc)
  end subroutine axissymlap_dlpn_blockmat_nmode_r64

end module axissymlap_specialquad_mod
