module axissym_physop_mod
  ! Physical-space SELF / OFF-DIAGONAL operators (mirror of test/omega/omega3_physical_op/AxiPhysMat.m
  ! and AxiOffDiagPhysMat.m).  Sits ABOVE the per-mode block builders: it calls the existing all-mode
  ! axissym{lap,stok}_*_blockmat_nmode_r64 routines, then assembles the modal stack into the same sparse
  ! factors MATLAB returns -- shipped across the mwrap boundary in CSR (ia=row_ptr, ja=col, a=val).
  use axistokes3d_mod, only: r64
  use axissymlap_specialquad_mod, only: axissymlap_slp_blockmat_nmode_r64,  axissymlap_slpn_blockmat_nmode_r64, &
                                        axissymlap_dlp_blockmat_nmode_r64,  axissymlap_dlpn_blockmat_nmode_r64
  use axissymstok_specialquad_mod, only: axissymstok_slp_blockmat_nmode_r64, axissymstok_slpn_blockmat_nmode_r64, &
                                         axissymstok_dlp_blockmat_nmode_r64, axissymstok_dlpn_blockmat_nmode_r64, &
                                         axissymstok_slppres_blockmat_nmode_r64, axissymstok_dlppres_blockmat_nmode_r64
  implicit none
  private
  public :: axissym_physmat_r64
  public :: axissym_offdiagphysmat_r64

contains

  ! ==================================================================================================
  ! AxiPhysMat: physical-space SELF operator, ALL modes in one call.
  !   iphys  1=laplace (nc=1) | 2=stokes (nc=3)        ilayer 1=slp 2=slpn 3=dlp 4=dlpn
  !   self => targets = source nodes (tx=sx, tnx=snx), nt = Nmer = p*np.
  ! Returns the five MATLAB factors [Ab, F, Finv, Abinv, T] in CSR (each N x N, N = rr*nphi):
  !   iAb/jAb/xAb : Ab     row_ptr(N+1) / col(nnzA) / val(nnzA)   nnzA = nphi*rr*rr   (rr entries/row)
  !   iAi/jAi/xAi : Abinv  (same shape as Ab)
  !   iF /jF /xF  : F      row_ptr(N+1) / col(nnzF) / val(nnzF)   nnzF = nphi*nphi*rr (nphi entries/row)
  !   iFi/jFi/xFi : Finv   (same shape as F)
  !   iT /jT /xT  : T      row_ptr(N+1) / col(nnzT) / val(nnzT)   nnzT = nc*nc*Nmer*nphi (nc entries/row)
  ! rr = nc*Nmer, Nmer = p*np, nphi = 2M+1, N = rr*nphi.  Indices 1-based.
  ! (NOTE this increment: iphys=1 laplace fully implemented; iphys=2 stokes pinv = next increment.)
  ! ==================================================================================================
  subroutine axissym_physmat_r64(iphys, ilayer, nc, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, &
                                 M, iside, iclosed, mu, R, xA, &
                                 iAb, jAb, xAb, iAi, jAi, xAi, iF, jF, xF, iFi, jFi, xFi, iT, jT, xT)
    integer(8),   intent(in)    :: iphys, ilayer, nc, p, np, M, iside, iclosed
    complex(r64), intent(in)    :: sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: sws(p*np), tpan(np+1), mu, R(3,3)
    ! xA = real(T*(Finv*Ab*F)*T') : assembled SELF operator, dense, LAB node-interleaved order (point-major)
    real(r64),    intent(inout) :: xA(nc*p*np*(2*M+1), nc*p*np*(2*M+1))
    ! CSR: i* = row_ptr (length N+1, N = nc*p*np*(2*M+1)) ; j* = col index ; x* = value
    real(r64),    intent(inout) :: iAb(nc*p*np*(2*M+1)+1), jAb((2*M+1)*(nc*p*np)**2)
    complex(r64), intent(inout) :: xAb((2*M+1)*(nc*p*np)**2)
    real(r64),    intent(inout) :: iAi(nc*p*np*(2*M+1)+1), jAi((2*M+1)*(nc*p*np)**2)
    complex(r64), intent(inout) :: xAi((2*M+1)*(nc*p*np)**2)
    real(r64),    intent(inout) :: iF(nc*p*np*(2*M+1)+1), jF((2*M+1)**2*(nc*p*np))
    complex(r64), intent(inout) :: xF((2*M+1)**2*(nc*p*np))
    real(r64),    intent(inout) :: iFi(nc*p*np*(2*M+1)+1), jFi((2*M+1)**2*(nc*p*np))
    complex(r64), intent(inout) :: xFi((2*M+1)**2*(nc*p*np))
    real(r64),    intent(inout) :: iT(nc*p*np*(2*M+1)+1), jT(nc*nc*p*np*(2*M+1))
    real(r64),    intent(inout) :: xT(nc*nc*p*np*(2*M+1))

    integer(8) :: Nmer, rr, nphi, Nful, md, n, k, i, jc, ka, kb, ii, aa, bb, iloc, node, gr, ptr, d
    real(r64)  :: twopi, phib, na
    complex(r64) :: ic
    real(r64), allocatable    :: Ar(:,:,:), RQall(:,:,:), tnphi0(:), Pnode(:,:,:), RQ(:,:,:)
    complex(r64), allocatable :: A(:,:,:), Ai(:,:,:), W(:,:), Winv(:,:), Cblk(:,:)
    real(r64) :: Q(3,3)

    Nmer = p*np; rr = nc*Nmer; nphi = 2*M+1; Nful = rr*nphi
    allocate(tnphi0(Nmer)); tnphi0 = 0.0_r64            ! SELF target: azimuthal normal n_theta = 0 (meridian normal)
    twopi = 2.0_r64*acos(-1.0_r64); ic = (0.0_r64, 1.0_r64)

    ! ---- per-mode block stack A(:,:,0..M) (complex; Laplace real promoted) ----
    allocate(A(rr, rr, M+1), Ai(rr, rr, M+1))
    if (iphys == 1) then
      allocate(Ar(Nmer, Nmer, M+1))
      select case (ilayer)
      case (1); call axissymlap_slp_blockmat_nmode_r64 (Nmer, sx,      p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, Ar)
      case (2); call axissymlap_slpn_blockmat_nmode_r64(Nmer, sx, snx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, Ar)
      case (3); call axissymlap_dlp_blockmat_nmode_r64 (Nmer, sx,      p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, Ar)
      case (4); call axissymlap_dlpn_blockmat_nmode_r64(Nmer, sx, snx, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, Ar)
      end select
      A = cmplx(Ar, 0.0_r64, r64)
      deallocate(Ar)
    else
      select case (ilayer)
      case (1); call axissymstok_slp_blockmat_nmode_r64 (Nmer, sx,      p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, mu, A)
      case (2); call axissymstok_slpn_blockmat_nmode_r64(Nmer, sx, snx, tnphi0, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, mu, A)
      case (3); call axissymstok_dlp_blockmat_nmode_r64 (Nmer, sx,      p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, mu, A)
      case (4); call axissymstok_dlpn_blockmat_nmode_r64(Nmer, sx, snx, tnphi0, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, mu, A)
      end select
    end if

    ! ---- per-mode inverse (Laplace: full-rank inv via Gauss-Jordan; Stokes: rank-deficient pinv via SVD) ----
    if (iphys == 1) then
      do md = 1, M+1
        call gj_inverse(rr, A(:,:,md), Ai(:,:,md))
      end do
    else
      do md = 1, M+1
        call pinv_svd(rr, A(:,:,md), Ai(:,:,md))
      end do
    end if

    ! ---- azimuthal DFT W, Winv  (W(ka,kb)=exp(-i n_ka phi_kb)/nphi ; Winv(kb,ka)=exp(i phi_kb n_ka)) ----
    allocate(W(nphi, nphi), Winv(nphi, nphi))
    do kb = 1, nphi
      phib = twopi*real(kb-1, r64)/real(nphi, r64)
      do ka = 1, nphi
        na = real(ka-1-M, r64)
        W(ka, kb)    = exp(-ic*na*phib)/real(nphi, r64)
        Winv(kb, ka) = real(nphi, r64)*conjg(W(ka, kb))          ! = exp(i phi_kb n_ka); no second exp
      end do
    end do

    ! ---- CSR: Ab, Abinv = blkdiag over n=-M..M (row gr=(k-1)*rr+ir has rr entries in block-cols) ----
    do gr = 1, Nful+1
      iAb(gr) = (gr-1)*rr + 1; iAi(gr) = (gr-1)*rr + 1
    end do
    ptr = 0
    do k = 1, nphi
      n = k - (M+1); md = abs(n) + 1
      do i = 1, rr                       ! row within block k
        do jc = 1, rr                    ! col within block k
          ptr = ptr + 1
          jAb(ptr) = (k-1)*rr + jc; jAi(ptr) = (k-1)*rr + jc
          if (n >= 0) then
            xAb(ptr) = A (i, jc, md); xAi(ptr) = Ai(i, jc, md)
          else
            xAb(ptr) = conjg(A (i, jc, md)); xAi(ptr) = conjg(Ai(i, jc, md))
          end if
        end do
      end do
    end do

    ! ---- CSR: F = kron(W, I_rr), Finv = kron(Winv, I_rr) (row gr=(ka-1)*rr+ir has nphi entries) ----
    do gr = 1, Nful+1
      iF(gr) = (gr-1)*nphi + 1; iFi(gr) = (gr-1)*nphi + 1
    end do
    ptr = 0
    do ka = 1, nphi
      do i = 1, rr
        do kb = 1, nphi
          ptr = ptr + 1
          jF(ptr)  = (kb-1)*rr + i; xF(ptr)  = W(ka, kb)
          jFi(ptr) = (kb-1)*rr + i; xFi(ptr) = Winv(ka, kb)
        end do
      end do
    end do

    ! ---- CSR: T = cyl->lab node-interleaved (nc=1 identity ; nc=3 per-azimuth R*Q fold) ----
    do gr = 1, Nful+1
      iT(gr) = (gr-1)*nc + 1
    end do
    if (nc == 1) then
      do gr = 1, Nful
        jT(gr) = gr; xT(gr) = 1.0_r64
      end do
    else
      allocate(RQall(3, 3, nphi))
      do ii = 1, nphi
        phib = twopi*real(ii-1, r64)/real(nphi, r64)
        Q = 0.0_r64
        Q(1,1) = cos(phib); Q(1,2) = -sin(phib)
        Q(2,1) = sin(phib); Q(2,2) =  cos(phib); Q(3,3) = 1.0_r64
        RQall(:,:,ii) = matmul(R, Q)
      end do
      ptr = 0
      do gr = 1, Nful                    ! lab node-interleaved row gr-1 = node*nc + bb
        bb   = mod(gr-1, nc)
        node = (gr-1)/nc
        ii   = node/Nmer + 1
        iloc = mod(node, Nmer)
        do aa = 0, nc-1
          ptr = ptr + 1
          jT(ptr) = (ii-1)*nc*Nmer + aa*Nmer + iloc + 1
          xT(ptr) = RQall(bb+1, aa+1, ii)
        end do
      end do
      deallocate(RQall)
    end if

    ! ---- dense assembled SELF operator xA = real(T*(Finv*Ab*F)*T') in LAB node-interleaved order ----
    !   azimuth-block-circulant P: Pcyl(d) = (1/nphi) sum_{n=-M..M} exp(i n d*dphi) A_n (A_n conj for n<0);
    !   reorder cyl-comp-major -> node-interleaved (Pnode); rotate per node by RQ = R*rotz(phi):
    !   xA[node-blk (i,I),(jc,J)] = RQ_i * Pnode(d=mod(i-jc,nphi))[I,J] * RQ_jc^T   (nc=1 -> RQ=1).
    allocate(Pnode(rr, rr, nphi), Cblk(rr, rr), RQ(nc, nc, nphi))
    do d = 0, nphi-1
      phib = twopi*real(d, r64)/real(nphi, r64)
      Cblk = (0.0_r64, 0.0_r64)
      do k = 1, nphi
        n = k - (M+1); md = abs(n) + 1
        if (n >= 0) then
          Cblk = Cblk + exp(ic*real(n, r64)*phib)*A(:,:,md)
        else
          Cblk = Cblk + exp(ic*real(n, r64)*phib)*conjg(A(:,:,md))
        end if
      end do
      do jc = 1, Nmer                    ! col meridian node (jloc)
        do bb = 0, nc-1                   ! col component
          do i = 1, Nmer                  ! row meridian node (iloc)
            do aa = 0, nc-1               ! row component
              Pnode((i-1)*nc+aa+1, (jc-1)*nc+bb+1, d+1) = &
                real(Cblk(aa*Nmer+i, bb*Nmer+jc), r64) / real(nphi, r64)
            end do
          end do
        end do
      end do
    end do
    do ii = 1, nphi                       ! per-azimuth rotation RQ = R*rotz(phi) (nc=1 -> 1)
      if (nc == 1) then
        RQ(1,1,ii) = 1.0_r64
      else
        phib = twopi*real(ii-1, r64)/real(nphi, r64)
        Q = 0.0_r64
        Q(1,1) = cos(phib); Q(1,2) = -sin(phib)
        Q(2,1) = sin(phib); Q(2,2) =  cos(phib); Q(3,3) = 1.0_r64
        RQ(:,:,ii) = matmul(R, Q)
      end if
    end do
    do i = 1, nphi                        ! azimuth row block
      do jc = 1, nphi                     ! azimuth col block
        d = mod(i - jc + nphi, nphi)
        do ii = 1, Nmer                   ! row meridian node I
          do k = 1, Nmer                  ! col meridian node J
            xA((i-1)*rr+(ii-1)*nc+1:(i-1)*rr+ii*nc, (jc-1)*rr+(k-1)*nc+1:(jc-1)*rr+k*nc) = &
              matmul(RQ(:,:,i), matmul(Pnode((ii-1)*nc+1:ii*nc, (k-1)*nc+1:k*nc, d+1), transpose(RQ(:,:,jc))))
          end do
        end do
      end do
    end do
    deallocate(Pnode, Cblk, RQ)

    deallocate(A, Ai, W, Winv)
  end subroutine axissym_physmat_r64

  ! ==================================================================================================
  ! AxiOffDiagPhysMat: physical-space OFF-DIAGONAL operator (one SOURCE particle -> ARBITRARY targets).
  !   Maps each target (tx=rho+iz, tphi=azimuth, frame Ct/Rt) into the SOURCE frame, classifies near/far
  !   (panel OR trapezoid-chi; iforce/=0 forces all near), and for NEAR targets builds the phase-baked
  !   block B = exp(i n th) .* A_n (A_n = the all-mode SLP value block at the mapped meridian txm).
  ! Out: B dense (nt x Nmer*nphi), near rows filled / far rows zero ; F = source DFT in CSR
  !      (ia=row_ptr(Nsrc+1)/ja/a, Nsrc=Nmer*nphi) ; near(nt) mask (1=near,0=far).
  !   operator action on a source physical density:  u_target(near) = real( B * (F * sigma) ).
  ! (NOTE this increment: only iphys=1,ilayer=1 laplace SLP implemented.)
  ! ==================================================================================================
  ! Full-3D LAB-geometry interface (mirror MATLAB axp_offdiagphysmat): TARGET given as Xt,Nt (LAB);
  ! source-frame reduction from (Cs,Rs).  SLPn (ilayer=2) uses ONE general-normal call -> sigma.n,
  ! INTERLOCKED B rows.  Other layers keep their old block-major internals (TODO).  sxlo/sxhi are
  ! derived in the mex .mw wrapper and passed in; the interlocked F reindex + At=real(B*F) are done there too.
  subroutine axissym_offdiagphysmat_r64(iphys, ilayer, nc, nt, ns, Xt, Ntg, Ct, Rt, &
                                        p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, Cs, Rs, &
                                        M, iside, iclosed, mu, iforce, &
                                        At, B, iF, jF, xF, near)
    integer(8),   intent(in)    :: iphys, ilayer, nc, nt, ns, p, np, M, iside, iclosed, iforce
    complex(r64), intent(in)    :: sx(p*np), snx(p*np), swxp(p*np), sxlo(np), sxhi(np)
    real(r64),    intent(in)    :: Xt(3,nt), Ntg(3,nt), Ct(3), Rt(3,3), Cs(3), Rs(3,3), sws(p*np), tpan(np+1), mu
    real(r64),    intent(inout) :: At(nc*nt, nc*ns)                          ! physical-space eval matrix, real(B*F), interlocked
    complex(r64), intent(inout) :: B(nc*nt, nc*p*np*(2*M+1))
    real(r64),    intent(inout) :: iF(nc*p*np*(2*M+1)+1), jF((2*M+1)**2*(nc*p*np))
    complex(r64), intent(inout) :: xF((2*M+1)**2*(nc*p*np))
    real(r64),    intent(inout) :: near(nt)

    integer(8) :: Nmer, ncN, nphi, i, j, k, n, kk, md, nnear, ka, kb, ptr, cnt, sp, comp, iblk, irow, iang
    integer(8), allocatable :: idn(:)
    real(r64)    :: twopi, gate, tol, chi_crit, panlen, rti, zti, chimin, chival, phib, na, thi, cc, ss
    real(r64)    :: nrho, nthe, Qs(3,3), vcyl(3), phis, cph, sph                      ! At source cyl->lab per-node rotation
    complex(r64) :: ic, e, ar, ap, az, vx, vy, vz, acc
    complex(r64) :: srr, str, szr, srz, stz, szz, stt, pp, Tr, Tt, Tz                 ! Stokes traction stress assembly
    real(r64)    :: xt3(3), nbf(3)
    real(r64),    allocatable :: th(:), thn(:), taz(:), nrhoa(:), nthea(:), nzza(:), Ar3(:,:,:), Av3(:,:,:)
    complex(r64), allocatable :: txm(:), txmn(:), tnxn(:), tnxn2(:), eone(:), eimg(:), W(:,:), Ast(:,:,:), Asz(:,:,:), Apr(:,:,:)
    real(r64),    allocatable :: tnp2(:)
    logical,      allocatable :: nearl(:)
    logical :: trac

    ! all 8 layers supported: Laplace SLP/SLPn/DLP/DLPn (nc=1), Stokes SLP/SLPn/DLP/DLPn (nc=3).
    if (iphys /= 1 .and. iphys /= 2) error stop 'axissym_offdiagphysmat_r64: bad iphys'
    if (ilayer < 1 .or. ilayer > 4)  error stop 'axissym_offdiagphysmat_r64: bad ilayer (1..4)'
    trac = (ilayer == 2 .or. ilayer == 4)         ! traction layers carry the target normal (tnx, tnphi)
    Nmer = p*np; ncN = nc*Nmer; nphi = 2*M+1
    twopi = 2.0_r64*acos(-1.0_r64); ic = (0.0_r64,1.0_r64)
    gate = 1.75_r64; tol = 1.0e-13_r64; chi_crit = cosh(log(1.0_r64/tol)/real(M, r64))

    ! ---- map all targets (LAB Xt) into the SOURCE frame: txm (meridian rho+iz) + th (azimuth) ----
    allocate(txm(nt), th(nt), nearl(nt))
    do i = 1, nt
      xt3 = matmul(transpose(Rs), Xt(:,i) - Cs)                                     ! target LAB -> source Cart
      txm(i) = cmplx(hypot(xt3(1), xt3(2)), xt3(3), r64)
      th(i)  = atan2(xt3(2), xt3(1))
    end do

    ! ---- classify near/far (panel rule OR trapezoid-chi rule); iforce/=0 => all near ----
    do i = 1, nt
      nearl(i) = (iforce /= 0)
    end do
    if (iforce == 0) then
      do k = 1, np                                  ! (1) PANEL rule
        panlen = 0.0_r64
        do j = 1, p; panlen = panlen + sws((k-1)*p+j); end do
        do i = 1, nt
          if (abs(txm(i)-sxlo(k)) + abs(txm(i)-sxhi(k)) < gate*panlen) nearl(i) = .true.
        end do
      end do
      do i = 1, nt                                  ! (2) TRAPEZOID chi rule
        if (.not. nearl(i)) then
          rti = real(txm(i)); zti = aimag(txm(i)); chimin = huge(1.0_r64)
          do j = 1, Nmer
            chival = 1.0_r64 + ((rti-real(sx(j)))**2 + (zti-aimag(sx(j)))**2)/(2.0_r64*rti*real(sx(j)))
            if (chival < chimin) chimin = chival
          end do
          if (chimin < chi_crit) nearl(i) = .true.
        end if
      end do
    end if
    do i = 1, nt
      if (nearl(i)) then; near(i) = 1.0_r64; else; near(i) = 0.0_r64; end if
    end do

    ! ---- gather near targets, call the all-mode block builder(s), scatter phase-baked into B ----
    nnear = count(nearl)
    B = (0.0_r64, 0.0_r64)
    if (nnear > 0) then
      allocate(idn(nnear), txmn(nnear), thn(nnear))
      cnt = 0
      do i = 1, nt
        if (nearl(i)) then; cnt = cnt+1; idn(cnt) = i; txmn(cnt) = txm(i); thn(cnt) = th(i); end if
      end do

      ! ---- traction layers: target normal (target frame) -> source-frame cyl components (nrho,nthe,nz) ----
      ! Laplace uses tnxn (meridian normal) + taz (n_theta/rho cross-term); Stokes contracts the full (nrho,nthe,nz).
      if (trac) then
        allocate(tnxn(nnear), taz(nnear), nrhoa(nnear), nthea(nnear), nzza(nnear)); cnt = 0
        do i = 1, nt
          if (.not. nearl(i)) cycle
          cnt = cnt + 1
          nbf = matmul(transpose(Rs), Ntg(:,i))                                    ! target LAB normal -> source frame (rotation only)
          thi = thn(cnt)
          nrho = nbf(1)*cos(thi) + nbf(2)*sin(thi)                                 ! meridian radial component
          nthe = -nbf(1)*sin(thi) + nbf(2)*cos(thi)                                ! azimuthal component
          tnxn(cnt) = cmplx(nrho, nbf(3), r64)                                     ! source-frame meridian normal (-> builder)
          taz(cnt)  = nthe / real(txmn(cnt))                                       ! azimuthal weight n_theta/rho (Laplace)
          nrhoa(cnt) = nrho; nthea(cnt) = nthe; nzza(cnt) = nbf(3)                  ! source-cyl normal (Stokes traction)
        end do
      end if

      if (iphys == 1) then                                    ! ===== LAPLACE (nc=1): scalar value (+ traction cross-term) =====
        allocate(Av3(nnear, Nmer, M+1))                       ! value block: SLP for slp/slpn, DLP for dlp/dlpn
        if (ilayer == 1 .or. ilayer == 2) then
          call axissymlap_slp_blockmat_nmode_r64(nnear, txmn, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, Av3)
        else
          call axissymlap_dlp_blockmat_nmode_r64(nnear, txmn, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, Av3)
        end if
        if (trac) then
          allocate(Ar3(nnear, Nmer, M+1))                     ! meridian-normal traction (gradient) block
          if (ilayer == 2) then
            call axissymlap_slpn_blockmat_nmode_r64(nnear, txmn, tnxn, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, Ar3)
          else
            call axissymlap_dlpn_blockmat_nmode_r64(nnear, txmn, tnxn, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, Ar3)
          end if
        end if
        do kk = 1, nphi
          n = kk - (M+1); md = abs(n) + 1
          do j = 1, Nmer
            do i = 1, nnear
              if (trac) then                                  ! traction = exp(i n th) (grad + i n (n_theta/rho) value)
                B(idn(i), (kk-1)*Nmer + j) = exp(ic*real(n,r64)*thn(i)) * &
                    (Ar3(i,j,md) + ic*real(n,r64)*taz(i)*Av3(i,j,md))
              else                                            ! value   = exp(i n th) A
                B(idn(i), (kk-1)*Nmer + j) = exp(ic*real(n,r64)*thn(i)) * Av3(i,j,md)
              end if
            end do
          end do
        end do
        deallocate(Av3); if (trac) deallocate(Ar3)
      else if (.not. trac) then                               ! ===== STOKES (nc=3) SLP/DLP value: cyl->src-Cart->lab(Rs) =====
        allocate(Ast(3*nnear, 3*Nmer, M+1))
        if (ilayer == 1) then
          call axissymstok_slp_blockmat_nmode_r64(nnear, txmn, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, mu, Ast)
        else
          call axissymstok_dlp_blockmat_nmode_r64(nnear, txmn, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, mu, Ast)
        end if
        do kk = 1, nphi
          n = kk - (M+1); md = abs(n) + 1
          do i = 1, nnear
            thi = thn(i); e = exp(ic*real(n,r64)*thi); cc = cos(thi); ss = sin(thi)
            do j = 1, 3*Nmer
              ar = Ast(i, j, md); ap = Ast(nnear+i, j, md); az = Ast(2*nnear+i, j, md)
              if (n < 0) then; ar = conjg(ar); ap = conjg(ap); az = conjg(az); end if
              vx = e*(cc*ar - ss*ap); vy = e*(ss*ar + cc*ap); vz = e*az
              B(nc*(idn(i)-1)+1, (kk-1)*ncN + j) = Rs(1,1)*vx + Rs(1,2)*vy + Rs(1,3)*vz   ! INTERLOCKED rows
              B(nc*(idn(i)-1)+2, (kk-1)*ncN + j) = Rs(2,1)*vx + Rs(2,2)*vy + Rs(2,3)*vz
              B(nc*(idn(i)-1)+3, (kk-1)*ncN + j) = Rs(3,1)*vx + Rs(3,2)*vy + Rs(3,3)*vz
            end do
          end do
        end do
        deallocate(Ast)
      else if (ilayer == 2) then                              ! ===== STOKES SLPn traction: ONE general-normal call -> sigma.n, INTERLOCKED B =====
        allocate(Ast(3*nnear,3*Nmer,M+1), tnxn2(nnear), tnp2(nnear))
        do i = 1, nnear
          tnxn2(i) = cmplx(nrhoa(i), nzza(i), r64); tnp2(i) = nthea(i)            ! target normal (source frame): meridian n_rho+i*n_z, azimuthal n_theta
        end do
        call axissymstok_slpn_blockmat_nmode_r64(nnear, txmn, tnxn2, tnp2, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, mu, Ast)
        do kk = 1, nphi
          n = kk - (M+1); md = abs(n) + 1
          do i = 1, nnear
            thi = thn(i); e = exp(ic*real(n,r64)*thi); cc = cos(thi); ss = sin(thi)
            do j = 1, 3*Nmer
              Tr = Ast(i,j,md); Tt = Ast(nnear+i,j,md); Tz = Ast(2*nnear+i,j,md)  ! t_rho, t_theta, t_z direct (sigma.n)
              if (n < 0) then; Tr = conjg(Tr); Tt = conjg(Tt); Tz = conjg(Tz); end if
              vx = e*(cc*Tr - ss*Tt); vy = e*(ss*Tr + cc*Tt); vz = e*Tz
              B(nc*(idn(i)-1)+1, (kk-1)*ncN + j) = Rs(1,1)*vx + Rs(1,2)*vy + Rs(1,3)*vz   ! INTERLOCKED rows
              B(nc*(idn(i)-1)+2, (kk-1)*ncN + j) = Rs(2,1)*vx + Rs(2,2)*vy + Rs(2,3)*vz
              B(nc*(idn(i)-1)+3, (kk-1)*ncN + j) = Rs(3,1)*vx + Rs(3,2)*vy + Rs(3,3)*vz
            end do
          end do
        end do
        deallocate(Ast, tnxn2, tnp2)
      else                                                    ! ===== STOKES DLPn traction: ONE general-normal call -> sigma.n, INTERLOCKED B =====
        allocate(Ast(3*nnear,3*Nmer,M+1), tnxn2(nnear), tnp2(nnear))
        do i = 1, nnear
          tnxn2(i) = cmplx(nrhoa(i), nzza(i), r64); tnp2(i) = nthea(i)            ! target normal (source frame): meridian n_rho+i*n_z, azimuthal n_theta
        end do
        call axissymstok_dlpn_blockmat_nmode_r64(nnear, txmn, tnxn2, tnp2, p, np, sx, snx, sws, swxp, tpan, sxlo, sxhi, M, iside, iclosed, mu, Ast)
        do kk = 1, nphi
          n = kk - (M+1); md = abs(n) + 1
          do i = 1, nnear
            thi = thn(i); e = exp(ic*real(n,r64)*thi); cc = cos(thi); ss = sin(thi)
            do j = 1, 3*Nmer
              Tr = Ast(i,j,md); Tt = Ast(nnear+i,j,md); Tz = Ast(2*nnear+i,j,md)  ! t_rho, t_theta, t_z direct (sigma.n)
              if (n < 0) then; Tr = conjg(Tr); Tt = conjg(Tt); Tz = conjg(Tz); end if
              vx = e*(cc*Tr - ss*Tt); vy = e*(ss*Tr + cc*Tt); vz = e*Tz
              B(nc*(idn(i)-1)+1, (kk-1)*ncN + j) = Rs(1,1)*vx + Rs(1,2)*vy + Rs(1,3)*vz   ! INTERLOCKED rows
              B(nc*(idn(i)-1)+2, (kk-1)*ncN + j) = Rs(2,1)*vx + Rs(2,2)*vy + Rs(2,3)*vz
              B(nc*(idn(i)-1)+3, (kk-1)*ncN + j) = Rs(3,1)*vx + Rs(3,2)*vy + Rs(3,3)*vz
            end do
          end do
        end do
        deallocate(Ast, tnxn2, tnp2)
      end if
      if (trac) deallocate(tnxn, taz, nrhoa, nthea, nzza)
      deallocate(idn, txmn, thn)
    end if

    ! ---- source DFT F = kron(W, I_Nmer) in CSR (Nsrc=Nmer*nphi rows, nphi entries/row) ----
    allocate(W(nphi, nphi))
    do kb = 1, nphi
      phib = twopi*real(kb-1, r64)/real(nphi, r64)
      do ka = 1, nphi
        na = real(ka-1-M, r64)
        W(ka, kb) = exp(-ic*na*phib)/real(nphi, r64)
      end do
    end do
    do i = 1, ncN*nphi+1
      iF(i) = (i-1)*nphi + 1
    end do
    ptr = 0
    do ka = 1, nphi
      do i = 1, ncN
        do kb = 1, nphi
          ptr = ptr + 1
          jF(ptr) = (kb-1)*ncN + i; xF(ptr) = W(ka, kb)
        end do
      end do
    end do
    ! ---- physical-space dense eval matrix At = real(B*F), source cols INTERLOCKED, then rotate source cyl -> LAB ----
    !   interlocked source col cc=(sp-1)*nc+comp; F picks B's angle-blocks: At_cyl(:,cc)=real(sum_ka B(:,(ka-1)*ncN+iblk) W(ka,iang)).
    !   for nc=3 the node's 3 cyl source comps are rotated to lab by Qs = Rs*rotz(phi_iang) (At = At_cyl*Qs', symmetric to the target Rs fold),
    !   so At is lab->lab: the caller feeds a LAB-Cartesian source density directly (no external per-node rotation).
    do sp = 1, ns
      iang = (sp-1)/Nmer + 1                                     ! source azimuthal block (1..nphi)
      j    = mod(sp-1, Nmer) + 1                                 ! source meridian node
      if (nc == 3) then
        phis = twopi*real(iang-1,r64)/real(nphi,r64); cph = cos(phis); sph = sin(phis)
        Qs(:,1) =  Rs(:,1)*cph + Rs(:,2)*sph                     ! Qs = Rs * rotz(phi) : source cyl -> lab
        Qs(:,2) = -Rs(:,1)*sph + Rs(:,2)*cph
        Qs(:,3) =  Rs(:,3)
      end if
      do irow = 1, nc*nt
        do comp = 1, nc
          iblk = (comp-1)*Nmer + j
          acc = (0.0_r64, 0.0_r64)
          do ka = 1, nphi
            acc = acc + B(irow, (ka-1)*ncN + iblk) * W(ka, iang)
          end do
          vcyl(comp) = real(acc)                                 ! source-cyl component
        end do
        if (nc == 3) then                                        ! At_lab(:,3k) = sum_m vcyl(m)*Qs(k,m)  (= At_cyl * Qs')
          At(irow, (sp-1)*3 + 1) = vcyl(1)*Qs(1,1) + vcyl(2)*Qs(1,2) + vcyl(3)*Qs(1,3)
          At(irow, (sp-1)*3 + 2) = vcyl(1)*Qs(2,1) + vcyl(2)*Qs(2,2) + vcyl(3)*Qs(2,3)
          At(irow, (sp-1)*3 + 3) = vcyl(1)*Qs(3,1) + vcyl(2)*Qs(3,2) + vcyl(3)*Qs(3,3)
        else
          At(irow, (sp-1)*nc + 1) = vcyl(1)                      ! nc=1 (Laplace): scalar source, no rotation
        end if
      end do
    end do
    deallocate(txm, th, nearl, W)
  end subroutine axissym_offdiagphysmat_r64

  ! --------------------------------------------------------------------------------------------------
  ! Dense complex inverse via Gauss-Jordan with partial pivoting (full-rank blocks only).
  ! --------------------------------------------------------------------------------------------------
  subroutine gj_inverse(n, Ain, Ainv)
    integer(8),   intent(in)  :: n
    complex(r64), intent(in)  :: Ain(n, n)
    complex(r64), intent(out) :: Ainv(n, n)
    complex(r64) :: aug(n, 2*n), piv, fac
    integer(8)   :: i, j, k, kmax
    real(r64)    :: amax
    aug = (0.0_r64, 0.0_r64)
    aug(:, 1:n) = Ain
    do i = 1, n
      aug(i, n+i) = (1.0_r64, 0.0_r64)
    end do
    do k = 1, n
      kmax = k; amax = abs(aug(k, k))
      do i = k+1, n
        if (abs(aug(i, k)) > amax) then; amax = abs(aug(i, k)); kmax = i; end if
      end do
      if (kmax /= k) then
        do j = 1, 2*n
          piv = aug(k, j); aug(k, j) = aug(kmax, j); aug(kmax, j) = piv
        end do
      end if
      piv = aug(k, k)
      aug(k, :) = aug(k, :)/piv
      do i = 1, n
        if (i /= k) then
          fac = aug(i, k)
          aug(i, :) = aug(i, :) - fac*aug(k, :)
        end if
      end do
    end do
    Ainv = aug(:, n+1:2*n)
  end subroutine gj_inverse

  ! --------------------------------------------------------------------------------------------------
  ! Dense complex Moore-Penrose pseudo-inverse via SVD (LAPACK zgesvd): rank-deficient safe (Stokes).
  ! Resolves against MATLAB's ILP64 libmwlapack at the MEX boundary, so the scalar args are plain
  ! integer(8) (= default under -fdefault-integer-8).  (OpenBLAS LAPACK is Flang-built and pulls an
  ! unavailable runtime symbol into the bundle, so it is NOT used here.)
  ! pinv = V S^+ U^H,  S^+ inverts the singular values > tol = n*eps*max(S) (MATLAB pinv tolerance).
  ! --------------------------------------------------------------------------------------------------
  subroutine pinv_svd(n, Ain, Ainv)
    integer(8),   intent(in)  :: n
    complex(r64), intent(in)  :: Ain(n, n)
    complex(r64), intent(out) :: Ainv(n, n)
    complex(r64), allocatable :: Acpy(:,:), U(:,:), VT(:,:), Vs(:,:), work(:)
    real(r64),    allocatable :: S(:), rwork(:)
    complex(r64) :: wq(1)
    real(r64)    :: tol
    integer(8)   :: lwork, info, i, k
    allocate(Acpy(n,n), U(n,n), VT(n,n), Vs(n,n), S(n), rwork(5*n))
    Acpy = Ain
    lwork = -1                                                   ! workspace query
    call zgesvd('A','A', n, n, Acpy, n, S, U, n, VT, n, wq, lwork, rwork, info)
    lwork = int(real(wq(1), r64), 8); allocate(work(lwork))
    call zgesvd('A','A', n, n, Acpy, n, S, U, n, VT, n, work, lwork, rwork, info)
    tol = real(n, r64)*epsilon(1.0_r64)*S(1)
    do k = 1, n                                                 ! Vs(:,k) = V(:,k)/s_k = conjg(VT(k,:))/s_k
      if (S(k) > tol) then
        do i = 1, n; Vs(i,k) = conjg(VT(k,i))/S(k); end do
      else
        do i = 1, n; Vs(i,k) = (0.0_r64, 0.0_r64); end do
      end if
    end do
    Ainv = matmul(Vs, conjg(transpose(U)))                      ! pinv = (V S^+) U^H
    deallocate(Acpy, U, VT, Vs, S, rwork, work)
  end subroutine pinv_svd

end module axissym_physop_mod
