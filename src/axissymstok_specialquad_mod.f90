module axissymstok_specialquad_mod
  ! ------------------------------------------------------------------
  ! Per-panel azimuthal-mode close-eval operator for the axisymmetric Stokes DLP / SLP
  ! velocity tensor.  Ports StokDLPAxiModalSpecialquad.m / StokAxiModalSpecialquad.m:
  ! assembles the per-mode interleaved 3nt x 3p block W by combining the kernel-split
  ! coefficients (axissymstok_kernelsplit_mod) with the Helsing singular quadrature
  ! (specialquad_mod) and the panel quadrature supplied by the caller (q->p Legendre
  ! projection P, upsampled weights gw).
  !
  ! Geometry is pre-evaluated by the host: the caller passes the upsampled source
  ! nodes y = Z(tt), tangents yp = Zp(tt), endpoints za=Z(a0), zb=Z(b0), and the
  ! parameter half-length hm = (b0-a0)/2 -- Fortran cannot take MATLAB handles.
  ! ------------------------------------------------------------------
  use axistokes3d_mod, only: r64
  use axissymstok_kernelsplit_mod, only: coef_dlp => axissymstok_kernelsplit_coef_dlp_r64, &
                                         coef_slp => axissymstok_kernelsplit_coef_slp_r64, &
                                         coef_slpn => axissymstok_kernelsplit_coef_slpn_r64
  use specialquad_mod, only: sdspecialquad_r64
  implicit none
  private
  public :: axissymstok_specialquad_dlp_r64
  public :: axissymstok_specialquad_slp_r64
  public :: axissymstok_specialquad_slpn_r64

contains

  ! ================================================================
  ! DLP close-eval block (port of StokDLPAxiModalSpecialquad).
  !   nt          number of targets
  !   z(nt)       complex targets
  !   p           output panel order; q = 2*p upsampled (be=2) source nodes
  !   hm          real, (b0-a0)/2  (parameter half-length of the panel)
  !   y(2p),yp(2p) complex source nodes and tangents (= Z(tt), Zp(tt))
  !   za, zb      complex panel endpoints (= Z(a0), Z(b0))
  !   gw(2p)      upsampled Gauss-Legendre weights on [-1,1]
  !   P(2p,p)     q->p Legendre projection matrix (= lege.pols(gx2,p-1).' * uk)
  !   M           max azimuthal mode
  !   mu          viscosity
  !   iside       1 = exterior ('e'), 0 = interior ('i')   (Sspecialquad/Dspecialquad branch)
  !   W(M+1,3nt,3p) REAL interleaved [rho;theta;zeta] block, per mode
  ! ================================================================
  subroutine axissymstok_specialquad_dlp_r64(nt, z, p, hm, y, yp, za, zb, gw, Pmat, M, mu, iside, W)
    integer(8),   intent(in)  :: nt, p, M, iside
    real(r64),    intent(in)  :: hm, mu, gw(2*p), Pmat(2*p,p)
    complex(r64), intent(in)  :: z(nt), y(2*p), yp(2*p), za, zb
    real(r64),    intent(inout) :: W(M+1,3*nt,3*p)

    integer(8) :: q, i, j, a, b, mm, row, col, kk
    real(r64)  :: pi, twopi, spdj
    complex(r64) :: ci, tauj, g3, g4, Azc
    real(r64)    :: g1, g2
    real(r64),    allocatable :: C1(:,:,:), C2(:,:,:)
    complex(r64), allocatable :: C3(:,:,:), C4(:,:,:), C5(:,:,:)
    real(r64),    allocatable :: As(:,:), A1(:,:), A2(:,:), A3(:,:), A4(:,:)
    complex(r64), allocatable :: Ad(:,:)
    real(r64),    allocatable :: arc(:), pref(:), Pk(:,:), blk(:,:), W2d(:,:)
    complex(r64), allocatable :: nv(:), wxp(:), tau(:)

    pi = acos(-1.0_r64); twopi = 2.0_r64*pi; ci = (0.0_r64,1.0_r64)
    q  = 2*p

    allocate(C1(M+1,3*nt,3*q), C2(M+1,3*nt,3*q))
    allocate(C3(M+1,3*nt,3*q), C4(M+1,3*nt,3*q), C5(M+1,3*nt,3*q))
    allocate(As(q,nt), Ad(q,nt), A1(q,nt), A2(q,nt), A3(q,nt), A4(q,nt))
    allocate(arc(q), pref(q), nv(q), wxp(q), tau(q))
    allocate(Pk(3*q,3*p), blk(3*nt,3*q), W2d(3*nt,3*p))

    ! geometry-derived per-node quantities (mirror the MATLAB panel setup)
    do j = 1, q
      spdj   = abs(yp(j))
      tau(j) = yp(j)/spdj
      nv(j)  = -ci*tau(j)
      arc(j) = spdj*gw(j)*hm
      wxp(j) = yp(j)*gw(j)*hm
      pref(j)= 3.0_r64*real(y(j),r64)/(4.0_r64*pi*mu)
    end do

    ! split coefficients (bare; pref applied below) and the singular quadrature
    call coef_dlp(nt, q, M, mu, real(z,r64), aimag(z), real(y,r64), aimag(y), &
                  real(nv,r64), aimag(nv), C1, C2, C3, C4, C5)
    call sdspecialquad_r64(nt, z, q, y, nv, wxp, za, zb, iside, As, Ad, A1, A2, A3, A4)

    ! q->p interleaved projection Pk = kron(P, I3)
    Pk = 0.0_r64
    do kk = 1, p
      do j = 1, q
        do b = 1, 3
          Pk(3*(j-1)+b, 3*(kk-1)+b) = Pmat(j,kk)
        end do
      end do
    end do

    ! per-mode assembly, then project q->p
    do mm = 0, M
      do i = 1, nt
        do a = 1, 3
          row = 3*(i-1)+a
          do j = 1, q
            tauj = tau(j)
            Azc  = cmplx(A1(j,i), -A2(j,i), r64)        ! Az = Re - i*(-Im) reconstructed
            do b = 1, 3
              col = 3*(j-1)+b
              g1 = pref(j)*C1(mm+1,row,col)
              g2 = pref(j)*C2(mm+1,row,col)
              g3 = pref(j)*C3(mm+1,row,col)
              g4 = pref(j)*C4(mm+1,row,col)
              blk(row,col) = twopi*As(j,i)*g1 + g2*arc(j) &
                   + real(twopi*ci*Ad(j,i)*(g3/tauj), r64) &
                   - real(twopi*ci*Azc*(g4/tauj), r64)
            end do
          end do
        end do
      end do
      W2d = matmul(blk, Pk)
      do b = 1, 3*p
        do row = 1, 3*nt
          W(mm+1,row,b) = W2d(row,b)
        end do
      end do
    end do

    deallocate(C1,C2,C3,C4,C5,As,Ad,A1,A2,A3,A4,arc,pref,nv,wxp,tau,Pk,blk,W2d)
  end subroutine axissymstok_specialquad_dlp_r64

  ! ================================================================
  ! SLP close-eval block (port of StokAxiModalSpecialquad).  Same panel inputs as the DLP.
  ! Coef map: coef_slp returns C1=log(G1), C2=smooth(G0), C3=Cauchy(Gc), each (M+1)x3nt x3q
  ! (fac=rho'/(4 mu) already folded).  Singular quad: As=Sspecialquad value (log),
  ! Ad=Dspecialquad value (Cauchy).  The Stokeslet Cauchy is handled by subtracting the
  ! canonical dyad Gcd (analytic, mode-indep) from Gc and adding the special-quad cauchyblk,
  ! both nonzero only in the 4 meridional entries (rho,rho),(rho,zeta),(zeta,rho),(zeta,zeta).
  !   blk = -2pi*G1.*KA + (Gc-Gcd+G0).*arc + cauchyblk ;  W = blk*Pk.
  ! ================================================================
  subroutine axissymstok_specialquad_slp_r64(nt, z, p, hm, y, yp, za, zb, gw, Pmat, M, mu, iside, W)
    integer(8),   intent(in)  :: nt, p, M, iside
    real(r64),    intent(in)  :: hm, mu, gw(2*p), Pmat(2*p,p)
    complex(r64), intent(in)  :: z(nt), y(2*p), yp(2*p), za, zb
    real(r64),    intent(inout) :: W(M+1,3*nt,3*p)

    integer(8) :: q, i, j, a, b, mm, row, col, kk
    real(r64)  :: pi, twopi, spdj, rti, ztimi, rsj, ysimj, w1, wz, w2a, cdij, g1, g0, gc, Gcdv, cauv, gs
    complex(r64) :: ci, Dmij, nxcj
    real(r64),    allocatable :: C1(:,:,:), C2(:,:,:), C3(:,:,:), C4(:,:,:), C5(:,:,:)
    real(r64),    allocatable :: As(:,:), A1(:,:), A2(:,:), A3(:,:), A4(:,:)
    complex(r64), allocatable :: Ad(:,:)
    real(r64),    allocatable :: arc(:), Pk(:,:), blk(:,:), W2d(:,:)
    complex(r64), allocatable :: nv(:), wxp(:)

    pi = acos(-1.0_r64); twopi = 2.0_r64*pi; ci = (0.0_r64,1.0_r64)
    q  = 2*p

    allocate(C1(M+1,3*nt,3*q), C2(M+1,3*nt,3*q))
    allocate(C3(M+1,3*nt,3*q), C4(M+1,3*nt,3*q), C5(M+1,3*nt,3*q))
    allocate(As(q,nt), Ad(q,nt), A1(q,nt), A2(q,nt), A3(q,nt), A4(q,nt))
    allocate(arc(q), nv(q), wxp(q))
    allocate(Pk(3*q,3*p), blk(3*nt,3*q), W2d(3*nt,3*p))

    do j = 1, q
      spdj   = abs(yp(j))
      nv(j)  = -ci*yp(j)/spdj
      arc(j) = spdj*gw(j)*hm
      wxp(j) = yp(j)*gw(j)*hm
    end do

    ! SLP split coefficients (fac folded) and the singular quadrature (As log, Ad Cauchy)
    call coef_slp(nt, q, M, mu, real(z,r64), aimag(z), real(y,r64), aimag(y), C1, C2, C3, C4, C5)
    call sdspecialquad_r64(nt, z, q, y, nv, wxp, za, zb, iside, As, Ad, A1, A2, A3, A4)

    Pk = 0.0_r64
    do kk = 1, p
      do j = 1, q
        do b = 1, 3
          Pk(3*(j-1)+b, 3*(kk-1)+b) = Pmat(j,kk)
        end do
      end do
    end do

    do mm = 0, M
      do i = 1, nt
        rti = real(z(i),r64); ztimi = aimag(z(i))
        do a = 1, 3
          row = 3*(i-1)+a
          do j = 1, q
            rsj = real(y(j),r64); ysimj = aimag(y(j))
            w1 = rsj - rti; wz = ysimj - ztimi; w2a = w1*w1 + wz*wz
            cdij = (rsj/(4.0_r64*mu)) / (pi*sqrt(rti*rsj))
            Dmij = Ad(j,i); nxcj = nv(j)
            do b = 1, 3
              col = 3*(j-1)+b
              g1 = C1(mm+1,row,col); g0 = C2(mm+1,row,col); gc = C3(mm+1,row,col)
              Gcdv = 0.0_r64; cauv = 0.0_r64
              if (a==1 .and. b==1) then
                Gcdv = cdij*w1*w1/w2a;  cauv = real(-twopi*Dmij*(w1/nxcj), r64)*cdij
              else if (a==1 .and. b==3) then
                Gcdv = cdij*w1*wz/w2a;  cauv = real(-twopi*Dmij*(wz/nxcj), r64)*cdij
              else if (a==3 .and. b==1) then
                Gcdv = cdij*w1*wz/w2a;  cauv = real(-twopi*Dmij*(ci*w1/nxcj), r64)*cdij
              else if (a==3 .and. b==3) then
                Gcdv = cdij*wz*wz/w2a;  cauv = real(-twopi*Dmij*(ci*wz/nxcj), r64)*cdij
              end if
              gs = gc - Gcdv + g0
              blk(row,col) = -twopi*g1*As(j,i) + gs*arc(j) + cauv
            end do
          end do
        end do
      end do
      W2d = matmul(blk, Pk)
      do b = 1, 3*p
        do row = 1, 3*nt
          W(mm+1,row,b) = W2d(row,b)
        end do
      end do
    end do

    deallocate(C1,C2,C3,C4,A1,A2,A3,A4,arc,nv,wxp,Pk,blk,W2d)
    deallocate(C5,As,Ad)
  end subroutine axissymstok_specialquad_slp_r64

  ! ================================================================
  ! SLPn close-eval block (port of StokSLPnAxiModalSpecialquad_v2).  TWO-BLOCK: W (meridional)
  ! and Wsw (swirl, unit n_theta).  Structurally identical to the DLP assembly, with:
  !   - the TARGET normal nx(nt) (per row) feeding coef_slpn (vs DLP's per-col source normal);
  !   - pref = -3 rho'/(4 pi) (mu-INDEPENDENT; DLP was +3 rho'/(4 pi mu));
  !   - W   : 2pi*As*g1 + g2*arc + Re(2pi i Ad g3/tau) - Re(2pi i Az g4/tau)   (mer, same as DLP)
  !   - Wsw : 2pi*As*g1 + g2*arc + Re(2pi i Ad g3/tau)                          (swirl, c4=0)
  ! The caller scales Wsw by n_theta and applies the FLIPPED parity mask.
  ! ================================================================
  subroutine axissymstok_specialquad_slpn_r64(nt, z, nx, p, hm, y, yp, za, zb, gw, Pmat, M, mu, iside, W, Wsw)
    integer(8),   intent(in)  :: nt, p, M, iside
    real(r64),    intent(in)  :: hm, mu, gw(2*p), Pmat(2*p,p)
    complex(r64), intent(in)  :: z(nt), nx(nt), y(2*p), yp(2*p), za, zb
    real(r64),    intent(inout) :: W(M+1,3*nt,3*p), Wsw(M+1,3*nt,3*p)

    integer(8) :: q, i, j, a, b, mm, row, col, kk
    real(r64)  :: pi, twopi, spdj
    complex(r64) :: ci, tauj, g3, Azc
    real(r64)    :: g1, g2
    real(r64),    allocatable :: C1(:,:,:), C2(:,:,:), C1s(:,:,:), C2s(:,:,:)
    complex(r64), allocatable :: C3(:,:,:), C4(:,:,:), C3s(:,:,:)
    real(r64),    allocatable :: As(:,:), A1(:,:), A2(:,:), A3(:,:), A4(:,:)
    complex(r64), allocatable :: Ad(:,:)
    real(r64),    allocatable :: arc(:), pref(:), Pk(:,:), blk(:,:), blks(:,:), W2d(:,:), Ws2d(:,:)
    complex(r64), allocatable :: nv(:), wxp(:), tau(:)

    pi = acos(-1.0_r64); twopi = 2.0_r64*pi; ci = (0.0_r64,1.0_r64)
    q  = 2*p

    allocate(C1(M+1,3*nt,3*q), C2(M+1,3*nt,3*q), C3(M+1,3*nt,3*q), C4(M+1,3*nt,3*q))
    allocate(C1s(M+1,3*nt,3*q), C2s(M+1,3*nt,3*q), C3s(M+1,3*nt,3*q))
    allocate(As(q,nt), Ad(q,nt), A1(q,nt), A2(q,nt), A3(q,nt), A4(q,nt))
    allocate(arc(q), pref(q), nv(q), wxp(q), tau(q))
    allocate(Pk(3*q,3*p), blk(3*nt,3*q), blks(3*nt,3*q), W2d(3*nt,3*p), Ws2d(3*nt,3*p))

    do j = 1, q
      spdj   = abs(yp(j))
      tau(j) = yp(j)/spdj
      nv(j)  = -ci*tau(j)
      arc(j) = spdj*gw(j)*hm
      wxp(j) = yp(j)*gw(j)*hm
      pref(j)= -3.0_r64*real(y(j),r64)/(4.0_r64*pi)
    end do

    ! split coefficients (TARGET normal nx; bare, pref applied below) + the singular quadrature
    call coef_slpn(nt, q, M, mu, real(z,r64), aimag(z), real(y,r64), aimag(y), &
                   real(nx,r64), aimag(nx), C1, C2, C3, C4, C1s, C2s, C3s)
    call sdspecialquad_r64(nt, z, q, y, nv, wxp, za, zb, iside, As, Ad, A1, A2, A3, A4)

    Pk = 0.0_r64
    do kk = 1, p
      do j = 1, q
        do b = 1, 3
          Pk(3*(j-1)+b, 3*(kk-1)+b) = Pmat(j,kk)
        end do
      end do
    end do

    do mm = 0, M
      do i = 1, nt
        do a = 1, 3
          row = 3*(i-1)+a
          do j = 1, q
            tauj = tau(j)
            Azc  = cmplx(A1(j,i), -A2(j,i), r64)
            do b = 1, 3
              col = 3*(j-1)+b
              ! meridional block W
              g1 = pref(j)*C1(mm+1,row,col)
              g2 = pref(j)*C2(mm+1,row,col)
              g3 = pref(j)*C3(mm+1,row,col)
              blk(row,col) = twopi*As(j,i)*g1 + g2*arc(j) &
                   + real(twopi*ci*Ad(j,i)*(g3/tauj), r64) &
                   - real(twopi*ci*Azc*(pref(j)*C4(mm+1,row,col)/tauj), r64)
              ! swirl block Wsw (unit n_theta; c4=0)
              g1 = pref(j)*C1s(mm+1,row,col)
              g2 = pref(j)*C2s(mm+1,row,col)
              g3 = pref(j)*C3s(mm+1,row,col)
              blks(row,col) = twopi*As(j,i)*g1 + g2*arc(j) &
                   + real(twopi*ci*Ad(j,i)*(g3/tauj), r64)
            end do
          end do
        end do
      end do
      W2d  = matmul(blk,  Pk)
      Ws2d = matmul(blks, Pk)
      do b = 1, 3*p
        do row = 1, 3*nt
          W(mm+1,row,b)   = W2d(row,b)
          Wsw(mm+1,row,b) = Ws2d(row,b)
        end do
      end do
    end do

    deallocate(C1,C2,C3,C4,C1s,C2s,C3s,As,Ad,A1,A2,A3,A4)
    deallocate(arc,pref,nv,wxp,tau,Pk,blk,blks,W2d,Ws2d)
  end subroutine axissymstok_specialquad_slpn_r64

end module axissymstok_specialquad_mod
