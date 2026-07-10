module axissym_physop_mod
  ! LEVEL 2 (toy/verification tier), PHYSICAL+MODAL setup masters.  Sits ABOVE the per-mode block
  ! builders (level 1): the masters take the full K-particle ragged geometry and return COMPLETE
  ! dense operators (no detection, every pair computed) -- axissym_physmat_setup assembles the
  ! physical-space operator, axissym_modemat_setup its Fourier-space twin (per-particle modal stacks).
  ! The physmat/offdiagphysmat workers below mirror test/omega/omega3_physical_op/AxiPhysMat.m and
  ! AxiOffDiagPhysMat.m; sparse factors ship across the mwrap boundary in CSR (ia=row_ptr, ja=col, a=val).
  use axistokes3d_mod, only: r64, gauss_r64, lagrange_interp_r64
  ! ---- LEVEL-1 WORKERS (kernel math): the all-mode block builders the masters dispatch to ----
  use axissymlap_specialquad_mod, only: &
       axissymlap_slp_blockmat_nmode_r64,   axissymlap_slpn_blockmat_nmode_r64,  &   ! Laplace SLP, SLPn
       axissymlap_dlp_blockmat_nmode_r64,   axissymlap_dlpn_blockmat_nmode_r64       ! Laplace DLP, DLPn
  use axissymstok_specialquad_mod, only: &
       axissymstok_slp_blockmat_nmode_r64,  axissymstok_slpn_blockmat_nmode_r64, &   ! Stokes  SLP, SLPn
       axissymstok_dlp_blockmat_nmode_r64,  axissymstok_dlpn_blockmat_nmode_r64, &   ! Stokes  DLP, DLPn
       axissymstok_slppres_blockmat_nmode_r64, axissymstok_dlppres_blockmat_nmode_r64 ! Stokes pressure
  implicit none
  private
  ! ---- MASTERS: the level-2 controlled surface ----
  public :: axissym_physmat_setup, axissym_modemat_setup
  ! ---- LEGACY WORKERS (kept byte-identical for their own mexes axp_physmat_mex /
  !      axp_offdiagphysmat_mex; physmat_setup no longer calls them) ----
  public :: axissym_physmat_r64, axissym_offdiagphysmat_r64
  ! ---- HELPERS (private): gj_inverse, pinv_svd, and the flag-free physmat_setup FOLD workers
  !      physmat_fold_self_r64 / physmat_fold_offdiag_r64 + phase-bake workers
  !      physmat_bake_lap_r64 / physmat_bake_lapn_r64 / physmat_bake_stok_r64 (verbatim copies
  !      of the fold stages of the legacy workers above) ----

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

  ! ==================================================================================================
  ! LEVEL-2 MASTER (physical half) -- APPROVED 2026-07-09, implemented.
  ! Mirror of axissym_close_setup MINUS gate/rball (L2 has no detection, every pair is computed),
  ! id -> (nrA, ncA, A): the COMPLETE dense operator (near close-eval + far analytic), toy tier.
  ! iinter: 1 self (K diagonal blocks, circulant) | 2 cross (off-diagonal blocks) | 3 eval.
  ! SAME dispatch shape as corr_setup/close_setup/modemat_setup: per particle the MASTER preps the
  ! targets (self: own meridian nodes; cross/eval: lab targets mapped into the source frame and
  ! near/far classified -- chi_crit trapezoid gate, iforce=0 semantics), then ONE nested
  ! select(ikernel) -> select(ilayer) dispatch to the PINNED *_blockmat_nmode_r64 workers builds
  ! the modal block stack, then the flag-free FOLD workers (physmat_bake_* + physmat_fold_self_r64
  ! / physmat_fold_offdiag_r64, verbatim copies of the fold stages of axissym_physmat_r64 /
  ! axissym_offdiagphysmat_r64, which stay untouched for their own mexes) assemble the dense block.
  ! A is real, LAB node-interleaved (per particle: azimuth-major 3D node order, nc comps/node).
  ! iinter=1: per-particle diagonal blocks nblk_k = nc*p_k*np_k*(2*pmodes(k)+1); off-diag untouched.
  ! iinter=2: targoff = per-particle OWN node ranges; each source kk fills the rows of every target
  !           OUTSIDE its own range (row of target g = nc*(g-1)) at its column block -- diagonal
  !           blocks untouched (caller combines with a self call, exactly like close_setup usage).
  ! iinter=3: arbitrary lab targets; every source sees ALL targets (targoff ignored); rows nc*ntarg.
  ! NOTE: A is NOT zeroed here -- the CALLER passes A pre-zeroed (the .mw wrapper allocates zeros),
  ! so the untouched blocks stay 0.  Rows of targets beyond the chi_crit trapezoid gate stay 0
  ! (modal far field, FMM territory).
  ! ==================================================================================================
  subroutine axissym_physmat_setup(ikernel, ilayer, params, iinter, K, p, np, pmodes, iside, iclosed, &
                                   geomoff, tpanoff, nsx, ntpan, sx, snx, sws, swxp, tpan, &
                                   ntarg, targoff, targ, targnx, Rm, Cc, nrA, ncA, A)
    integer(8),   intent(in)    :: ikernel, ilayer, iinter, K, iside, iclosed, nsx, ntpan, ntarg, nrA, ncA
    integer(8),   intent(in)    :: p(K), np(K), pmodes(K), geomoff(K+1), tpanoff(K+1), targoff(K+1)
    complex(r64), intent(in)    :: params, sx(nsx), snx(nsx), swxp(nsx)
    real(r64),    intent(in)    :: sws(nsx), tpan(ntpan), targ(3,ntarg), targnx(3,ntarg), Rm(3,3,K), Cc(3,K)
    real(r64),    intent(inout) :: A(nrA, ncA)
    integer(8) :: nc, kk, pk, npk, pmk, Nk, nangk, nblk, go, tpo, co, ns3, ntk, ntt, nnear, i, j, g, cnt, kp
    real(r64)  :: mu, rle(2), gate, tol, chi_crit, panlen, rti, zti, chimin, chival, thi, nrho, nthe
    real(r64)  :: xt3(3), nbf(3)
    logical    :: trac
    complex(r64) :: ec2(2)
    real(r64),    allocatable :: tglp(:), wglp(:), Dp(:,:), Lend(:,:)
    complex(r64), allocatable :: sxlo_k(:), sxhi_k(:)
    ! modal block stack (pinned blockmat_nmode output) + per-iinter target prep
    real(r64),    allocatable :: Akr(:,:,:), Akv(:,:,:), tnphik(:), th(:), thn(:), taz(:)
    complex(r64), allocatable :: Akz(:,:,:), Ac(:,:,:), txk(:), tnxk(:), txm(:), Bk(:,:)
    integer(8),   allocatable :: tsel(:), idn(:)
    real(r64),    allocatable :: Xtk(:,:), Ngk(:,:), Atk(:,:), xAk(:,:)
    logical,      allocatable :: nearl(:)
    select case (ikernel)                                        ! 1 lap (nc=1) / 2 stok (nc=3)
    case (1_8); nc = 1; mu = 0.0_r64
    case (2_8); nc = 3; mu = real(params, r64)
    case default; print *, 'axissym_physmat_setup: unknown ikernel', ikernel; stop
    end select
    if (ilayer == 5_8 .or. ilayer == 6_8) then
      print *, 'axissym_physmat_setup: ilayer', ilayer, ' pressure not wired in physmat_setup'; stop
    end if
    if (ilayer < 1_8 .or. ilayer > 4_8) then; print *, 'axissym_physmat_setup: bad ilayer', ilayer; stop; end if
    if (iinter < 1_8 .or. iinter > 3_8) then; print *, 'axissym_physmat_setup: bad iinter', iinter; stop; end if
    ! ---- fail loud on dimension mismatch: cols = nc*(total 3D source nodes); rows by iinter ----
    ns3 = 0
    do kk = 1, K
      ns3 = ns3 + p(kk)*np(kk)*(2*pmodes(kk)+1)
    end do
    if (ncA /= nc*ns3) then; print *, 'axissym_physmat_setup: ncA', ncA, ' /= nc*sum(Nk*nang)', nc*ns3; stop; end if
    if (iinter == 1_8) then
      if (nrA /= nc*ns3) then; print *, 'axissym_physmat_setup: self nrA', nrA, ' /= nc*sum(Nk*nang)', nc*ns3; stop; end if
    else
      if (nrA /= nc*ntarg) then; print *, 'axissym_physmat_setup: nrA', nrA, ' /= nc*ntarg', nc*ntarg; stop; end if
    end if
    if (iinter == 2_8) then                                      ! cross: targets ARE the surface nodes
      do kk = 1, K
        if (targoff(kk+1) - targoff(kk) /= p(kk)*np(kk)*(2*pmodes(kk)+1)) then
          print *, 'axissym_physmat_setup: cross targoff block', kk, ' size', targoff(kk+1)-targoff(kk), &
                   ' /= particle 3D node count', p(kk)*np(kk)*(2*pmodes(kk)+1); stop
        end if
      end do
    end if
    rle(1) = -1.0_r64; rle(2) = 1.0_r64
    trac = (ilayer == 2_8 .or. ilayer == 4_8)                    ! traction layers carry the target normal
    gate = 1.75_r64; tol = 1.0e-13_r64                           ! offdiag near/far gates (verbatim offdiagphysmat)
    co = 0
    do kk = 1, K                                                 ! ragged extraction, corr_setup-style
      pk = p(kk); npk = np(kk); pmk = pmodes(kk); Nk = pk*npk; nangk = 2*pmk+1
      nblk = nc*Nk*nangk; go = geomoff(kk); tpo = tpanoff(kk)
      ! ---- panel endpoints sxlo_k/sxhi_k from the generating curve sx (Lagrange endpoint interp) ----
      allocate(tglp(pk), wglp(pk), Dp(pk,pk), Lend(2,pk), sxlo_k(npk), sxhi_k(npk))
      call gauss_r64(pk, tglp, wglp, Dp)
      call lagrange_interp_r64(pk, tglp, 2_8, rle, Lend)
      do j = 1, npk
        ec2 = matmul(Lend, sx(go+(j-1)*pk : go+j*pk-1))
        sxlo_k(j) = ec2(1); sxhi_k(j) = ec2(2)
      end do
      ! ---- TARGET PREP (select iinter): SELF = own meridian nodes ; CROSS/EVAL = lab targets ----
      ! mapped into the SOURCE frame + near/far classified (chi_crit trapezoid gate, iforce=0),
      ! near targets gathered (idn/txk/thn) with traction normals -> source-frame cyl components.
      ntk = 0; ntt = 0; nnear = 0
      select case (iinter)
      case (1_8)                                                 ! SELF: tnphi = 0 (meridian normal)
        ntt = Nk
        allocate(txk(Nk), tnxk(Nk), tnphik(Nk))
        txk = sx(go:go+Nk-1); tnxk = snx(go:go+Nk-1); tnphik = 0.0_r64
      case (2_8, 3_8)                                            ! CROSS (non-own targets) / EVAL (all targets)
        if (iinter == 2_8) then
          ntk = ntarg - (targoff(kk+1) - targoff(kk))
        else
          ntk = ntarg
        end if
        if (ntk > 0) then
          allocate(tsel(ntk), Xtk(3,ntk), Ngk(3,ntk))
          cnt = 0
          do g = 1, ntarg
            if (iinter == 2_8 .and. g >= targoff(kk) .and. g < targoff(kk+1)) cycle
            cnt = cnt + 1; tsel(cnt) = g; Xtk(:,cnt) = targ(:,g); Ngk(:,cnt) = targnx(:,g)
          end do
          ! ---- map all targets (LAB Xtk) into the SOURCE frame: txm (meridian rho+iz) + th (azimuth) ----
          allocate(txm(ntk), th(ntk), nearl(ntk))
          do i = 1, ntk
            xt3 = matmul(transpose(Rm(:,:,kk)), Xtk(:,i) - Cc(:,kk))               ! target LAB -> source Cart
            txm(i) = cmplx(hypot(xt3(1), xt3(2)), xt3(3), r64)
            th(i)  = atan2(xt3(2), xt3(1))
          end do
          ! ---- classify near/far (panel rule OR trapezoid-chi rule); iforce=0 semantics ----
          chi_crit = cosh(log(1.0_r64/tol)/real(pmk, r64))
          do i = 1, ntk
            nearl(i) = .false.
          end do
          do kp = 1, npk                                ! (1) PANEL rule
            panlen = 0.0_r64
            do j = 1, pk; panlen = panlen + sws(go-1+(kp-1)*pk+j); end do
            do i = 1, ntk
              if (abs(txm(i)-sxlo_k(kp)) + abs(txm(i)-sxhi_k(kp)) < gate*panlen) nearl(i) = .true.
            end do
          end do
          do i = 1, ntk                                 ! (2) TRAPEZOID chi rule
            if (.not. nearl(i)) then
              rti = real(txm(i)); zti = aimag(txm(i)); chimin = huge(1.0_r64)
              do j = 1, Nk
                chival = 1.0_r64 + ((rti-real(sx(go-1+j)))**2 + (zti-aimag(sx(go-1+j)))**2)/(2.0_r64*rti*real(sx(go-1+j)))
                if (chival < chimin) chimin = chival
              end do
              if (chimin < chi_crit) nearl(i) = .true.
            end if
          end do
          ! ---- gather near targets (+ traction normals: LAB -> source-frame cyl components) ----
          nnear = count(nearl)
          allocate(Bk(nc*ntk, nblk)); Bk = (0.0_r64, 0.0_r64)
          if (nnear > 0) then
            allocate(idn(nnear), txk(nnear), thn(nnear))
            cnt = 0
            do i = 1, ntk
              if (nearl(i)) then; cnt = cnt+1; idn(cnt) = i; txk(cnt) = txm(i); thn(cnt) = th(i); end if
            end do
            if (trac) then
              allocate(tnxk(nnear), taz(nnear), tnphik(nnear)); cnt = 0
              do i = 1, ntk
                if (.not. nearl(i)) cycle
                cnt = cnt + 1
                nbf = matmul(transpose(Rm(:,:,kk)), Ngk(:,i))                      ! target LAB normal -> source frame (rotation only)
                thi = thn(cnt)
                nrho = nbf(1)*cos(thi) + nbf(2)*sin(thi)                           ! meridian radial component
                nthe = -nbf(1)*sin(thi) + nbf(2)*cos(thi)                          ! azimuthal component
                tnxk(cnt)   = cmplx(nrho, nbf(3), r64)                             ! source-frame meridian normal (-> builder)
                taz(cnt)    = nthe / real(txk(cnt))                                ! azimuthal weight n_theta/rho (Laplace)
                tnphik(cnt) = nthe                                                 ! azimuthal normal n_theta (Stokes)
              end do
            end if
          end if
          ntt = nnear
        end if
      end select
      ! ---- MODAL BLOCK STACK: ONE nested (ikernel, ilayer) dispatch to the PINNED blockmat workers ----
      ! (same shape as modemat_setup); cross/eval branches also phase-bake the stack into Bk.
      if (iinter == 1_8 .or. nnear > 0) then
        if (ikernel == 1_8) then                                 ! LAPLACE: real Akr(ntt, Nk, pmk+1)
          allocate(Akr(ntt, Nk, pmk+1))
          select case (ilayer)
          case (1_8)                                             ! ---- SLP ----
            call axissymlap_slp_blockmat_nmode_r64 (ntt, txk,       pk, npk, sx(go), snx(go), sws(go), swxp(go), &
                                                    tpan(tpo), sxlo_k, sxhi_k, pmk, iside, iclosed, Akr)
            if (iinter /= 1_8) call physmat_bake_lap_r64(ntk, ntt, Nk, pmk, idn, thn, Akr, Bk)
          case (2_8)                                             ! ---- SLPn (+ SLP value block for the n_theta/rho cross-term) ----
            call axissymlap_slpn_blockmat_nmode_r64(ntt, txk, tnxk, pk, npk, sx(go), snx(go), sws(go), swxp(go), &
                                                    tpan(tpo), sxlo_k, sxhi_k, pmk, iside, iclosed, Akr)
            if (iinter /= 1_8) then
              allocate(Akv(ntt, Nk, pmk+1))
              call axissymlap_slp_blockmat_nmode_r64(ntt, txk,      pk, npk, sx(go), snx(go), sws(go), swxp(go), &
                                                     tpan(tpo), sxlo_k, sxhi_k, pmk, iside, iclosed, Akv)
              call physmat_bake_lapn_r64(ntk, ntt, Nk, pmk, idn, thn, taz, Akr, Akv, Bk)
              deallocate(Akv)
            end if
          case (3_8)                                             ! ---- DLP ----
            call axissymlap_dlp_blockmat_nmode_r64 (ntt, txk,       pk, npk, sx(go), snx(go), sws(go), swxp(go), &
                                                    tpan(tpo), sxlo_k, sxhi_k, pmk, iside, iclosed, Akr)
            if (iinter /= 1_8) call physmat_bake_lap_r64(ntk, ntt, Nk, pmk, idn, thn, Akr, Bk)
          case (4_8)                                             ! ---- DLPn (+ DLP value block for the n_theta/rho cross-term) ----
            call axissymlap_dlpn_blockmat_nmode_r64(ntt, txk, tnxk, pk, npk, sx(go), snx(go), sws(go), swxp(go), &
                                                    tpan(tpo), sxlo_k, sxhi_k, pmk, iside, iclosed, Akr)
            if (iinter /= 1_8) then
              allocate(Akv(ntt, Nk, pmk+1))
              call axissymlap_dlp_blockmat_nmode_r64(ntt, txk,      pk, npk, sx(go), snx(go), sws(go), swxp(go), &
                                                     tpan(tpo), sxlo_k, sxhi_k, pmk, iside, iclosed, Akv)
              call physmat_bake_lapn_r64(ntk, ntt, Nk, pmk, idn, thn, taz, Akr, Akv, Bk)
              deallocate(Akv)
            end if
          case default; print *, 'axissym_physmat_setup: bad ilayer', ilayer; stop
          end select
        else                                                     ! STOKES: complex Akz(3*ntt, 3*Nk, pmk+1)
          allocate(Akz(3*ntt, 3*Nk, pmk+1))
          select case (ilayer)
          case (1_8); call axissymstok_slp_blockmat_nmode_r64 (ntt, txk,                pk, npk, sx(go), snx(go), sws(go), swxp(go), &
                                                               tpan(tpo), sxlo_k, sxhi_k, pmk, iside, iclosed, mu, Akz)
          case (2_8); call axissymstok_slpn_blockmat_nmode_r64(ntt, txk, tnxk, tnphik, pk, npk, sx(go), snx(go), sws(go), swxp(go), &
                                                               tpan(tpo), sxlo_k, sxhi_k, pmk, iside, iclosed, mu, Akz)
          case (3_8); call axissymstok_dlp_blockmat_nmode_r64 (ntt, txk,                pk, npk, sx(go), snx(go), sws(go), swxp(go), &
                                                               tpan(tpo), sxlo_k, sxhi_k, pmk, iside, iclosed, mu, Akz)
          case (4_8); call axissymstok_dlpn_blockmat_nmode_r64(ntt, txk, tnxk, tnphik, pk, npk, sx(go), snx(go), sws(go), swxp(go), &
                                                               tpan(tpo), sxlo_k, sxhi_k, pmk, iside, iclosed, mu, Akz)
          case default; print *, 'axissym_physmat_setup: bad ilayer', ilayer; stop
          end select
          if (iinter /= 1_8) call physmat_bake_stok_r64(ntk, ntt, Nk, pmk, idn, thn, Rm(:,:,kk), Akz, Bk)
        end if
      end if
      ! ---- FOLD (flag-free workers): modal stack -> dense physical block, place into A ----
      select case (iinter)
      case (1_8)                                                 ! SELF: circulant fold -> diagonal block
        allocate(Ac(nc*Nk, nc*Nk, pmk+1), xAk(nblk, nblk))
        if (ikernel == 1_8) then
          Ac = cmplx(Akr, 0.0_r64, r64)
        else
          Ac = Akz
        end if
        call physmat_fold_self_r64(nc, Nk, pmk, Rm(:,:,kk), Ac, xAk)
        A(co+1:co+nblk, co+1:co+nblk) = xAk
        deallocate(Ac, xAk)
      case (2_8, 3_8)                                            ! CROSS/EVAL: DFT + rotation fold -> target rows
        if (ntk > 0) then
          allocate(Atk(nc*ntk, nblk))
          call physmat_fold_offdiag_r64(nc, ntk, Nk*nangk, Nk, pmk, Rm(:,:,kk), Bk, Atk)
          do i = 1, ntk                                          ! place rows at the target's global row offset
            g = tsel(i)
            A(nc*(g-1)+1:nc*g, co+1:co+nblk) = Atk(nc*(i-1)+1:nc*i, :)
          end do
          deallocate(Atk)
        end if
      end select
      co = co + nblk
      if (allocated(Akr))    deallocate(Akr)
      if (allocated(Akz))    deallocate(Akz)
      if (allocated(txk))    deallocate(txk)
      if (allocated(tnxk))   deallocate(tnxk)
      if (allocated(tnphik)) deallocate(tnphik)
      if (allocated(taz))    deallocate(taz)
      if (allocated(idn))    deallocate(idn)
      if (allocated(thn))    deallocate(thn)
      if (allocated(txm))    deallocate(txm)
      if (allocated(th))     deallocate(th)
      if (allocated(nearl))  deallocate(nearl)
      if (allocated(tsel))   deallocate(tsel)
      if (allocated(Xtk))    deallocate(Xtk)
      if (allocated(Ngk))    deallocate(Ngk)
      if (allocated(Bk))     deallocate(Bk)
      deallocate(tglp, wglp, Dp, Lend, sxlo_k, sxhi_k)
    end do
    if (co /= ncA) then; print *, 'axissym_physmat_setup: col total', co, ' /= ncA', ncA; stop; end if
  end subroutine axissym_physmat_setup

  ! ==================================================================================================
  ! physmat_setup FOLD workers (private, flag-free): VERBATIM COPIES of the fold/assembly stages of
  ! axissym_physmat_r64 / axissym_offdiagphysmat_r64 (the legacy workers stay untouched, their mexes
  ! axp_physmat_mex / axp_offdiagphysmat_mex depend on them byte-identical).  Each takes the prebuilt
  ! modal block stack + (nc, geometry, frames) only -- no kernel/layer flags.
  ! ==================================================================================================

  ! ---- SELF fold: per-mode block stack A(:,:,0..M) -> dense assembled SELF operator xA ----
  ! (verbatim: the "dense assembled SELF operator" stage of axissym_physmat_r64)
  subroutine physmat_fold_self_r64(nc, Nmer, M, R, A, xA)
    integer(8),   intent(in)    :: nc, Nmer, M
    real(r64),    intent(in)    :: R(3,3)
    complex(r64), intent(in)    :: A(nc*Nmer, nc*Nmer, M+1)
    real(r64),    intent(inout) :: xA(nc*Nmer*(2*M+1), nc*Nmer*(2*M+1))
    integer(8) :: rr, nphi, md, n, k, i, jc, ii, aa, bb, d
    real(r64)  :: twopi, phib
    complex(r64) :: ic
    real(r64)  :: Q(3,3)
    real(r64), allocatable    :: Pnode(:,:,:), RQ(:,:,:)
    complex(r64), allocatable :: Cblk(:,:)
    rr = nc*Nmer; nphi = 2*M+1
    twopi = 2.0_r64*acos(-1.0_r64); ic = (0.0_r64, 1.0_r64)
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
  end subroutine physmat_fold_self_r64

  ! ---- OFFDIAG phase-bake, LAPLACE value layers (slp/dlp): B(near row) = exp(i n th) A ----
  ! (verbatim: the non-traction Laplace scatter loop of axissym_offdiagphysmat_r64)
  subroutine physmat_bake_lap_r64(nt, nnear, Nmer, M, idn, thn, Av3, B)
    integer(8),   intent(in)    :: nt, nnear, Nmer, M
    integer(8),   intent(in)    :: idn(nnear)
    real(r64),    intent(in)    :: thn(nnear), Av3(nnear, Nmer, M+1)
    complex(r64), intent(inout) :: B(nt, Nmer*(2*M+1))
    integer(8) :: nphi, kk, n, md, i, j
    complex(r64) :: ic
    nphi = 2*M+1; ic = (0.0_r64, 1.0_r64)
    do kk = 1, nphi
      n = kk - (M+1); md = abs(n) + 1
      do j = 1, Nmer
        do i = 1, nnear
          B(idn(i), (kk-1)*Nmer + j) = exp(ic*real(n,r64)*thn(i)) * Av3(i,j,md)
        end do
      end do
    end do
  end subroutine physmat_bake_lap_r64

  ! ---- OFFDIAG phase-bake, LAPLACE traction layers (slpn/dlpn): grad block + i n (n_theta/rho) value ----
  ! (verbatim: the traction Laplace scatter loop of axissym_offdiagphysmat_r64)
  subroutine physmat_bake_lapn_r64(nt, nnear, Nmer, M, idn, thn, taz, Ar3, Av3, B)
    integer(8),   intent(in)    :: nt, nnear, Nmer, M
    integer(8),   intent(in)    :: idn(nnear)
    real(r64),    intent(in)    :: thn(nnear), taz(nnear), Ar3(nnear, Nmer, M+1), Av3(nnear, Nmer, M+1)
    complex(r64), intent(inout) :: B(nt, Nmer*(2*M+1))
    integer(8) :: nphi, kk, n, md, i, j
    complex(r64) :: ic
    nphi = 2*M+1; ic = (0.0_r64, 1.0_r64)
    do kk = 1, nphi
      n = kk - (M+1); md = abs(n) + 1
      do j = 1, Nmer
        do i = 1, nnear                                 ! traction = exp(i n th) (grad + i n (n_theta/rho) value)
          B(idn(i), (kk-1)*Nmer + j) = exp(ic*real(n,r64)*thn(i)) * &
              (Ar3(i,j,md) + ic*real(n,r64)*taz(i)*Av3(i,j,md))
        end do
      end do
    end do
  end subroutine physmat_bake_lapn_r64

  ! ---- OFFDIAG phase-bake, STOKES (all four layers): cyl -> src-Cart -> lab(Rs), INTERLOCKED rows ----
  ! (verbatim: the Stokes scatter loop of axissym_offdiagphysmat_r64 -- identical for slp/slpn/dlp/dlpn)
  subroutine physmat_bake_stok_r64(nt, nnear, Nmer, M, idn, thn, Rs, Ast, B)
    integer(8),   intent(in)    :: nt, nnear, Nmer, M
    integer(8),   intent(in)    :: idn(nnear)
    real(r64),    intent(in)    :: thn(nnear), Rs(3,3)
    complex(r64), intent(in)    :: Ast(3*nnear, 3*Nmer, M+1)
    complex(r64), intent(inout) :: B(3*nt, 3*Nmer*(2*M+1))
    integer(8) :: ncN, nphi, kk, n, md, i, j
    real(r64)  :: thi, cc, ss
    complex(r64) :: ic, e, ar, ap, az, vx, vy, vz
    ncN = 3*Nmer; nphi = 2*M+1; ic = (0.0_r64, 1.0_r64)
    do kk = 1, nphi
      n = kk - (M+1); md = abs(n) + 1
      do i = 1, nnear
        thi = thn(i); e = exp(ic*real(n,r64)*thi); cc = cos(thi); ss = sin(thi)
        do j = 1, 3*Nmer
          ar = Ast(i, j, md); ap = Ast(nnear+i, j, md); az = Ast(2*nnear+i, j, md)
          if (n < 0) then; ar = conjg(ar); ap = conjg(ap); az = conjg(az); end if
          vx = e*(cc*ar - ss*ap); vy = e*(ss*ar + cc*ap); vz = e*az
          B(3*(idn(i)-1)+1, (kk-1)*ncN + j) = Rs(1,1)*vx + Rs(1,2)*vy + Rs(1,3)*vz   ! INTERLOCKED rows
          B(3*(idn(i)-1)+2, (kk-1)*ncN + j) = Rs(2,1)*vx + Rs(2,2)*vy + Rs(2,3)*vz
          B(3*(idn(i)-1)+3, (kk-1)*ncN + j) = Rs(3,1)*vx + Rs(3,2)*vy + Rs(3,3)*vz
        end do
      end do
    end do
  end subroutine physmat_bake_stok_r64

  ! ---- OFFDIAG fold: phase-baked B -> dense physical eval matrix At = real(B*F), lab->lab ----
  ! (verbatim: the source-DFT W + At stage of axissym_offdiagphysmat_r64, minus the CSR emission)
  subroutine physmat_fold_offdiag_r64(nc, nt, ns, Nmer, M, Rs, B, At)
    integer(8),   intent(in)    :: nc, nt, ns, Nmer, M
    real(r64),    intent(in)    :: Rs(3,3)
    complex(r64), intent(in)    :: B(nc*nt, nc*Nmer*(2*M+1))
    real(r64),    intent(inout) :: At(nc*nt, nc*ns)
    integer(8) :: ncN, nphi, j, ka, kb, sp, comp, iblk, irow, iang
    real(r64)  :: twopi, phib, na, phis, cph, sph, Qs(3,3), vcyl(3)
    complex(r64) :: ic, acc
    complex(r64), allocatable :: W(:,:)
    ncN = nc*Nmer; nphi = 2*M+1
    twopi = 2.0_r64*acos(-1.0_r64); ic = (0.0_r64,1.0_r64)
    allocate(W(nphi, nphi))
    do kb = 1, nphi
      phib = twopi*real(kb-1, r64)/real(nphi, r64)
      do ka = 1, nphi
        na = real(ka-1-M, r64)
        W(ka, kb) = exp(-ic*na*phib)/real(nphi, r64)
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
    deallocate(W)
  end subroutine physmat_fold_offdiag_r64

  ! ==================================================================================================
  ! LEVEL-2 MASTER (modal half) -- Fourier-space twin of axissym_physmat_setup: IDENTICAL signature;
  ! only two differences:
  !   (1) iinter=2 (cross) FAILS LOUD -- inter-particle coupling mixes azimuthal modes;
  !   (2) A carries the mode dimension: A(nrA, ncA, maxval(pmodes)+1), complex for BOTH kernels.
  ! iinter=1: K decoupled per-particle SELF modal stacks (mode-by-mode solves).
  ! iinter=3: modal eval blocks, lab targets mapped into each particle meridian (Rm^T(x-Cc)).
  ! Dispatches to the *_blockmat_nmode_r64 workers.
  ! ilayer 1=slp 2=slpn 3=dlp 4=dlpn | 5=slppres 6=dlppres (STOKES ONLY, ikernel=2): pressure rows
  ! are SCALAR -- per particle the block is (ntk) x (3*Nk), so nrA = sum(ntk) while ncA stays 3*nsx.
  ! NOTE: A is NOT zeroed here (assumed-size 3rd dim) -- the CALLER passes A pre-zeroed (the .mw
  ! wrapper allocates complex(zeros(...))); higher modes of low-pmode particles stay 0.
  ! ==================================================================================================
  subroutine axissym_modemat_setup(ikernel, ilayer, params, iinter, K, p, np, pmodes, iside, iclosed, &
                                   geomoff, tpanoff, nsx, ntpan, sx, snx, sws, swxp, tpan, &
                                   ntarg, targoff, targ, targnx, Rm, Cc, nrA, ncA, A)
    integer(8),   intent(in)    :: ikernel, ilayer, iinter, K, iside, iclosed, nsx, ntpan, ntarg, nrA, ncA
    integer(8),   intent(in)    :: p(K), np(K), pmodes(K), geomoff(K+1), tpanoff(K+1), targoff(K+1)
    complex(r64), intent(in)    :: params, sx(nsx), snx(nsx), swxp(nsx)
    real(r64),    intent(in)    :: sws(nsx), tpan(ntpan), targ(3,ntarg), targnx(3,ntarg), Rm(3,3,K), Cc(3,K)
    complex(r64), intent(inout) :: A(nrA, ncA, *)                ! 3rd dim = maxval(pmodes)+1
    integer(8) :: nc, ncr, kk, pk, npk, pmk, Nk, go, tpo, ntk, ntsum, ro, co, i, j, m, tg
    real(r64)  :: mu, th, rle(2), v(3), nv(3)
    real(r64),    allocatable :: tglp(:), wglp(:), Dp(:,:), Lend(:,:), tnphik(:), Akr(:,:,:)
    complex(r64), allocatable :: sxlo_k(:), sxhi_k(:), txk(:), tnxk(:), Akz(:,:,:)
    complex(r64) :: ec2(2)
    select case (ikernel)
    case (1_8); nc = 1; mu = 0.0_r64
    case (2_8); nc = 3; mu = real(params, r64)
    case default; print *, 'axissym_modemat_setup: unknown ikernel', ikernel; stop
    end select
    if (ilayer < 1_8 .or. ilayer > 6_8) then; print *, 'axissym_modemat_setup: bad ilayer', ilayer; stop; end if
    if (ilayer >= 5_8 .and. ikernel /= 2_8) then
      print *, 'axissym_modemat_setup: ilayer', ilayer, ' (pressure) requires ikernel=2 (Stokes)'; stop
    end if
    if (iinter == 2_8) then; print *, 'axissym_modemat_setup: cross couples modes, no Fourier meaning'; stop; end if
    if (iinter /= 1_8 .and. iinter /= 3_8) then; print *, 'axissym_modemat_setup: bad iinter', iinter; stop; end if
    ! ---- fail loud on dimension mismatch: rows = ncr*(total targets), cols = nc*(total source nodes) ----
    ! ncr = per-target row count: nc components, EXCEPT ilayer 5/6 pressure -> SCALAR rows (ncr=1)
    ncr = nc; if (ilayer >= 5_8) ncr = 1_8
    if (iinter == 1_8) then; ntsum = nsx; else; ntsum = targoff(K+1) - targoff(1); end if
    if (nrA /= ncr*ntsum) then; print *, 'axissym_modemat_setup: nrA', nrA, ' /= ncr*ntsum', ncr*ntsum; stop; end if
    if (ncA /= nc*nsx)    then; print *, 'axissym_modemat_setup: ncA', ncA, ' /= nc*nsx', nc*nsx; stop; end if
    rle(1) = -1.0_r64; rle(2) = 1.0_r64
    ro = 0
    do kk = 1, K                                                 ! ragged extraction, corr_setup-style
      pk = p(kk); npk = np(kk); pmk = pmodes(kk); Nk = pk*npk; go = geomoff(kk); tpo = tpanoff(kk)
      ! ---- panel endpoints sxlo_k/sxhi_k from the generating curve sx (Lagrange endpoint interp) ----
      allocate(tglp(pk), wglp(pk), Dp(pk,pk), Lend(2,pk), sxlo_k(npk), sxhi_k(npk))
      call gauss_r64(pk, tglp, wglp, Dp)
      call lagrange_interp_r64(pk, tglp, 2_8, rle, Lend)
      do j = 1, npk
        ec2 = matmul(Lend, sx(go+(j-1)*pk : go+j*pk-1))
        sxlo_k(j) = ec2(1); sxhi_k(j) = ec2(2)
      end do
      ! ---- targets: SELF = own meridian nodes ; EVAL = lab targets mapped into this particle frame ----
      select case (iinter)
      case (1_8)                                                 ! SELF: tnphi = 0 (meridian normal)
        ntk = Nk
        allocate(txk(ntk), tnxk(ntk), tnphik(ntk))
        txk = sx(go:go+Nk-1); tnxk = snx(go:go+Nk-1); tnphik = 0.0_r64
      case (3_8)                                                 ! EVAL: v = Rm^T (x - Cc)
        ntk = targoff(kk+1) - targoff(kk)
        allocate(txk(ntk), tnxk(ntk), tnphik(ntk))
        do i = 1, ntk
          tg = targoff(kk) - 1 + i
          v  = matmul(transpose(Rm(:,:,kk)), targ(:,tg) - Cc(:,kk))
          txk(i) = cmplx(hypot(v(1), v(2)), v(3), r64)
          th = atan2(v(2), v(1))                                 ! target azimuth in the particle frame
          nv = matmul(transpose(Rm(:,:,kk)), targnx(:,tg))       ! lab normal -> particle frame (rotation only)
          tnxk(i)   = cmplx(nv(1)*cos(th) + nv(2)*sin(th), nv(3), r64)   ! meridian normal n_rho + i n_z
          tnphik(i) = -nv(1)*sin(th) + nv(2)*cos(th)             ! azimuthal normal n_theta
        end do
      end select
      ! ---- dispatch (ikernel, ilayer) -> blockmat_nmode worker; place at (ro, co) ----
      co = nc*(geomoff(kk) - 1)
      if (ikernel == 1_8) then                                   ! LAPLACE: real Ak(ntk, Nk, pmk+1)
        allocate(Akr(ntk, Nk, pmk+1))
        select case (ilayer)
        case (1_8); call axissymlap_slp_blockmat_nmode_r64 (ntk, txk,       pk, npk, sx(go), snx(go), sws(go), swxp(go), &
                                                            tpan(tpo), sxlo_k, sxhi_k, pmk, iside, iclosed, Akr)
        case (2_8); call axissymlap_slpn_blockmat_nmode_r64(ntk, txk, tnxk, pk, npk, sx(go), snx(go), sws(go), swxp(go), &
                                                            tpan(tpo), sxlo_k, sxhi_k, pmk, iside, iclosed, Akr)
        case (3_8); call axissymlap_dlp_blockmat_nmode_r64 (ntk, txk,       pk, npk, sx(go), snx(go), sws(go), swxp(go), &
                                                            tpan(tpo), sxlo_k, sxhi_k, pmk, iside, iclosed, Akr)
        case (4_8); call axissymlap_dlpn_blockmat_nmode_r64(ntk, txk, tnxk, pk, npk, sx(go), snx(go), sws(go), swxp(go), &
                                                            tpan(tpo), sxlo_k, sxhi_k, pmk, iside, iclosed, Akr)
        case default; print *, 'axissym_modemat_setup: bad ilayer', ilayer; stop
        end select
        do m = 1, pmk+1
          A(ro+1:ro+ntk, co+1:co+Nk, m) = cmplx(Akr(:,:,m), 0.0_r64, r64)
        end do
        deallocate(Akr)
      else                                                       ! STOKES: complex Ak(ncr*ntk, 3*Nk, pmk+1)
        allocate(Akz(ncr*ntk, 3*Nk, pmk+1))                      ! ncr=3 velocity rows | ncr=1 pressure (scalar) rows
        select case (ilayer)
        case (1_8); call axissymstok_slp_blockmat_nmode_r64 (ntk, txk,                pk, npk, sx(go), snx(go), sws(go), swxp(go), &
                                                             tpan(tpo), sxlo_k, sxhi_k, pmk, iside, iclosed, mu, Akz)
        case (2_8); call axissymstok_slpn_blockmat_nmode_r64(ntk, txk, tnxk, tnphik, pk, npk, sx(go), snx(go), sws(go), swxp(go), &
                                                             tpan(tpo), sxlo_k, sxhi_k, pmk, iside, iclosed, mu, Akz)
        case (3_8); call axissymstok_dlp_blockmat_nmode_r64 (ntk, txk,                pk, npk, sx(go), snx(go), sws(go), swxp(go), &
                                                             tpan(tpo), sxlo_k, sxhi_k, pmk, iside, iclosed, mu, Akz)
        case (4_8); call axissymstok_dlpn_blockmat_nmode_r64(ntk, txk, tnxk, tnphik, pk, npk, sx(go), snx(go), sws(go), swxp(go), &
                                                             tpan(tpo), sxlo_k, sxhi_k, pmk, iside, iclosed, mu, Akz)
        ! ---- 5/6: SLP/DLP PRESSURE (mu-free, no tnx/tnphi): scalar rows Ak(ntk, 3*Nk, pmk+1) ----
        case (5_8); call axissymstok_slppres_blockmat_nmode_r64(ntk, txk,            pk, npk, sx(go), snx(go), sws(go), swxp(go), &
                                                             tpan(tpo), sxlo_k, sxhi_k, pmk, iside, iclosed, Akz)
        case (6_8); call axissymstok_dlppres_blockmat_nmode_r64(ntk, txk,            pk, npk, sx(go), snx(go), sws(go), swxp(go), &
                                                             tpan(tpo), sxlo_k, sxhi_k, pmk, iside, iclosed, Akz)
        case default; print *, 'axissym_modemat_setup: bad ilayer', ilayer; stop
        end select
        do m = 1, pmk+1
          A(ro+1:ro+ncr*ntk, co+1:co+3*Nk, m) = Akz(:,:,m)
        end do
        deallocate(Akz)
      end if
      ro = ro + ncr*ntk
      deallocate(tglp, wglp, Dp, Lend, sxlo_k, sxhi_k, txk, tnxk, tnphik)
    end do
    if (ro /= nrA) then; print *, 'axissym_modemat_setup: row total', ro, ' /= nrA', nrA; stop; end if
  end subroutine axissym_modemat_setup

end module axissym_physop_mod
