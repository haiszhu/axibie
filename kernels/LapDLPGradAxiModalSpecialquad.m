function W = LapDLPGradAxiModalSpecialquad(z, ntn, a0, b0, p, gx, gw, Z, Zp, M, iself)
% Per-panel azimuthal-mode close-eval block for the axisymmetric Laplace DLPn D'
% (hypersingular: target- and source-normal double derivative of G).
%   W = LapDLPGradAxiModalSpecialquad(z, ntn, a0, b0, p, gx, gw, Z, Zp, M, iself) -> (M+1) x p
% ntn = target meridian unit normal [n_rho; n_z].  FOUR-PART split: log (Sspecialquad) +
% smooth (Gauss) + Cauchy (Dspecialquad value) + deriv-Cauchy (Dspecialquad Az).  Per-mode
% coefficients (G1,c2,c3,c4) from the validated util axissymlap_dlpn_kernel_split4 (chunkie
% convention), rescaled to physical (1/(2pi rho_src)).  c3,c4 already carry the target/source
% normals.  iself = self-node index, or [] for near-but-not-self.

hm=(b0-a0)/2; mid=(a0+b0)/2; tt=mid+hm*gx;
y=Z(tt); yp=Zp(tt); spd=abs(yp).'; tau=yp./abs(yp); nv=-1i*tau;     % source tangent/normal
rp=real(y).';
src=struct('r',[real(y(:))';imag(y(:))'],'n',[real(nv(:))';imag(nv(:))']);
tgt=struct('r',[real(z);imag(z)],'n',ntn(:));

ss.x=y(:); ss.wxp=yp(:).*gw(:)*hm; ss.nx=nv(:); tg.x=z;
Aspq=Sspecialquad(tg,ss,Z(a0),Z(b0),'e'); Aspq=Aspq(:).';
[Acau,Azcau]=Dspecialquad(tg,ss,Z(a0),Z(b0),'e'); Acau=Acau(:).'; Azcau=Azcau(:).';
if ~isempty(iself)                                                % self node: principal value
    [Aci,Azi]=Dspecialquad(tg,ss,Z(a0),Z(b0),'i'); Acau=(Acau+Aci(:).')/2; Azcau=(Azcau+Azi(:).')/2;
end
taur=tau(:).'; gwr=gw.'; convp=1./(2*pi*rp);                        % chunkie -> physical, per node

W=zeros(M+1,p);
for m=0:M
    [G1,c2,c3,c4]=axissymlap_dlpn_kernel_split4(src, tgt, m);      % chunkie conv, 1 x p
    G1p=G1.*convp; c2p=c2.*convp; c3p=c3.*convp; c4p=c4.*convp;
    logblk     = -2*pi*Aspq.*G1p;
    smoothblk  =  c2p.*(spd.*gwr*hm);
    cauchyblk  =  real( 2*pi*1i*Acau .* (c3p./taur) );
    dcauchyblk = -real( 2*pi*1i*Azcau .* (c4p./taur) );           % Az = d/dx[1/(x-y)] = -1/(x-y)^2
    W(m+1,:)   = logblk + smoothblk + cauchyblk + dcauchyblk;
end
end

function [G1, c2, c3, c4] = axissymlap_dlpn_kernel_split4(srcinfo, targinfo, m)
% ALL-MODE FOUR-PART split of the per-mode axisym Laplace DLPn (hypersingular) kernel:
%
%   K = log(rhat)*G1 + c2 + Re[c3/(r-r')] + Re[c4/(r-r')^2],   (r-r') = dr + i dz.
%
% G1,c2 REAL; c3,c4 COMPLEX; all SMOOTH at coincidence (the explicit 1/(r-r') and
% 1/(r-r')^2 carry the Cauchy and hypersingular singularities, log(rhat) the log).
%
% Derived from the validated all-mode six-piece (axissymlap_dlpn_kernel_split) by:
%  (i)  (n_x.d)(n_y.d)/rhat^4 = 1/2 Re[n nu/(r-r')^2] + 1/2 (n.nu)/rhat^2; the real
%       Hadamard (n.nu)/rhat^2 cancels (H_nn + 1/2 H_dd = 0)  =>  c4 = 1/2 H_dd n_x n_y;
%  (ii) the HIDDEN Cauchy mm'/rhat^2 = Re[m_t nu/(r-r')] (m_t=n_x.d, m'=n_y.d) buried in
%       the six-piece smooth block G0 (its 2 F' Phi/(rho rho' rhat^2) term) is pulled out
%       into c3; what remains, c2 = G0 - Re[hidden/(r-r')], is genuinely smooth (cf. the
%       Stokes split's G0 = K - log G1 - Gc by subtraction).
% This is the all-mode generalization of the 0th-mode elliptic recipe in
% laplace_dprime_axi.tex / LapDGradAxiSpecialMat1pan (function DGradCoeffs).

[Hnn,Hdd,Cr,Cz,G1,G0] = axissymlap_dlpn_kernel_split(srcinfo, targinfo, m);

src_r=srcinfo.r(1,:); src_z=srcinfo.r(2,:); nrs=srcinfo.n(1,:); nzs=srcinfo.n(2,:);
targ_r=targinfo.r(1,:); targ_z=targinfo.r(2,:); nrt=targinfo.n(1,:).'; nzt=targinfo.n(2,:).';
ns=numel(src_r); nt=numel(targ_r);
rt=repmat(targ_r.',1,ns); zt=repmat(targ_z.',1,ns);
rs=repmat(src_r,nt,1);    zs=repmat(src_z,nt,1);
nx1=repmat(nrt,1,ns); nx2=repmat(nzt,1,ns); ny1=repmat(nrs,nt,1); ny2=repmat(nzs,nt,1);
dr=rt-rs; dz=zt-zs; rhat2=dr.^2+dz.^2; rr=rt.*rs; B=1./sqrt(rr);
x=-rhat2./(4*rr);                                            % 2F1 argument = (1-chi)/2
Fp=0.5*(m^2-0.25).*hyp2f1_series(1.5-m,1.5+m,2,x,80);       % F' = P'_{m-1/2}

mt=nx1.*dr+nx2.*dz;  ms=ny1.*dr+ny2.*dz;                    % n_x.d , n_y.d
nt_c=nx1+1i*nx2; nu_c=ny1+1i*ny2;                          % complex meridian normals

hid = (rs/pi).*(B.*Fp./rr).*mt;                            % hidden-Cauchy coeff (smooth, O(rhat))
c3  = (Cr + 1i*Cz) + hid.*nu_c;                            % Cauchy
c4  = 0.5*Hdd.*(nt_c.*nu_c);                               % deriv-Cauchy

hidval = (rs/pi).*(B.*Fp./rr).*(mt.*ms./rhat2);           % = Re[hid*nu/(r-r')]
hidval(rhat2==0)=0;
c2  = G0 - hidval;                                          % smooth remainder
end

% ===========================================================================
function s = hyp2f1_series(a,b,c,x,N)
s = ones(size(x)); term = ones(size(x));
for n = 0:N-1
    term = term .* ((a+n).*(b+n)./((c+n)*(n+1))) .* x;
    s = s + term;
end
end

function [H_nn, H_dd, Cr, Cz, G1, G0] = axissymlap_dlpn_kernel_split(srcinfo, targinfo, m)
% Split the per-mode axisymmetric Laplace DLPn (hypersingular) kernel into
%
%   K = (n_x.n_y)/rhat^2 * H_nn + (n_x.d)(n_y.d)/rhat^4 * H_dd
%       + (rho-rho')/rhat^2 * Cr + (z-z')/rhat^2 * Cz + log(rhat) * G1 + G0,
%
% d = (rho-rho', z-z'), n_x = target meridian normal, n_y = source meridian normal.
% Mirrors axissymlap_slpn_kernel_split: the SAME Q,Q',Q'' values are regrouped so the
% explicit factors carry the 1/rhat^2 (hyper), 1/rhat (Cauchy) and log singularities
% while H_nn, H_dd, Cr, Cz, G0 stay bounded at coincidence.  Matches
% axissymlap_dlpn_kernel's convention (K = +2*pi*rho_src * physical D').
%
% Method: with u = chi-1, F = P_{m-1/2}(chi), R the bounded remainder in
%   Q = -1/2 log(u) F + R,
%   Q'  = -1/2 F/u - 1/2 log(u) F' + R',
%   Q'' =  1/2 F/u^2 - F'/u - 1/2 log(u) F'' + R'',
% the Brick-1 Hessian contraction is linear in (Q,Q',Q''), so
%   K = K_log + K_rat + K_R,
%   K_log = -1/2 log(u) K_F  (K_F built from F,F',F''),
%   K_rat built from (0, -1/2 F/u, 1/2 F/u^2 - F'/u),
%   K_R   built from (R,R',R'').
% log(u) = 2 log(rhat) - log(2 rho rho') gives G1 = -K_F and +1/2 log(2 rho rho') K_F
% into G0; K_R is bounded -> G0; the singular content of K_rat is pulled into the
% explicit factors with bounded coefficients (see below).

src_r  = srcinfo.r(1,:);   src_z  = srcinfo.r(2,:);
nrs = srcinfo.n(1,:);      nzs = srcinfo.n(2,:);
targ_r = targinfo.r(1,:);  targ_z = targinfo.r(2,:);
nrt = targinfo.n(1,:).';   nzt = targinfo.n(2,:).';

ns = numel(src_r);  nt = numel(targ_r);
rt = repmat(targ_r.',1,ns);  zt = repmat(targ_z.',1,ns);
rs = repmat(src_r,   nt,1);  zs = repmat(src_z,   nt,1);
NRT = repmat(nrt,1,ns); NZT = repmat(nzt,1,ns);
NRS = repmat(nrs,nt,1); NZS = repmat(nzs,nt,1);

dr = rt - rs;  dz = zt - zs;  rhat2 = dr.^2 + dz.^2;
u = rhat2 ./ (2*rt.*rs);                       % chi - 1

[F, Fp, Fpp]  = F21_and_2derivatives(u, m);
[Rr, Rp, Rpp] = remainder_and_2derivatives(u, m);

% --- log part: K_log = -1/2 log(u) K_F ------------------------------------
KF = contract_QQdQdd(rt,zt,rs,zs,NRT,NZT,NRS,NZS, F, Fp, Fpp);
G1 = -KF;
G0_log = 0.5*log(2*rt.*rs) .* KF;

% --- remainder part: K_R bounded -> G0 ------------------------------------
KR = contract_QQdQdd(rt,zt,rs,zs,NRT,NZT,NRS,NZS, Rr, Rp, Rpp);

% --- singular factors (all closed form, cancellation-free) -----------------
% From Q'' = 1/2 F/u^2 + ... and n_x.grad chi = a/(rho rho')+O(rhat^2),
% n_y.grad' chi = -b/(rho rho')+O(rhat^2), with 1/u^2 = 4 rho^2 rho'^2/rhat^4:
B = 1./sqrt(rt.*rs);
G = (rs/(2*pi)) .* B .* F;            % common factor
H_dd = -2*G;                          % coeff of (n_x.d)(n_y.d)/rhat^4  (= -(rs/pi) B F)
H_nn =  G;                            % coeff of (n_x.n_y)/rhat^2 (coincidence value)
Cr = G .* NRT.*NRS .* dr ./ (2*rt.*rs);          % coeff of (rho-rho')/rhat^2
Cz = G .* ( NZT.*NRS./(2*rs) - NRT.*NZS./(2*rt) );% coeff of (z-z')/rhat^2

% --- bounded remainder of the rational core (CLOSED FORM, no 1/u^2) --------
% Px.*Py./rhat2 is well-conditioned (Px,Py = O(rhat)), so no cancellation.
aa = NRT.*dr + NZT.*dz;  bb = NRS.*dr + NZS.*dz;            % n_x.d, n_y.d
Px =  aa./(rt.*rs) - NRT.*rhat2./(2*rt.^2.*rs);            % n_x.grad chi
Py = -bb./(rt.*rs) - NRS.*rhat2./(2*rt.*rs.^2);            % n_y.grad' chi
G0_rat = (rs/(2*pi)) .* B .* ( -F.*NRT.*NRS./(2*rt.*rs) ...
                               - 2*Fp.*rt.*rs.*(Px.*Py./rhat2) );

G0 = G0_log + KR + G0_rat;
end

% ===========================================================================
function K = contract_QQdQdd(rt,zt,rs,zs,NRT,NZT,NRS,NZS, Q,Qd,Qdd)
% K = (rs/2pi) n_x^T Hess(B Q) n_y, the Brick-1 DLPn contraction
% (axissymlap_dlpn_kernel) but with Q,Qd,Qdd passed in so callers can substitute the
% F-parts / R-parts of the toroidal harmonic.
dz = zt - zs;  rhat2 = (rt-rs).^2 + dz.^2;
B  = 1./sqrt(rt.*rs);  Br = -B./(2*rt);  Brp = -B./(2*rs);  Brrp = B./(4*rt.*rs);
chir  =  (rt-rs)./(rt.*rs) - rhat2./(2*rt.^2.*rs);
chiz  =  dz./(rt.*rs);
chirp = -(rt-rs)./(rt.*rs) - rhat2./(2*rt.*rs.^2);
chizp = -dz./(rt.*rs);
chirrp = -1./rs.^2 + (rt-rs)./(rt.^2.*rs) + rhat2./(2*rt.^2.*rs.^2);
chirzp =  dz./(rt.^2.*rs);
chizrp = -dz./(rt.*rs.^2);
chizzp = -1./(rt.*rs);
Grrp = Brrp.*Q + Br.*Qd.*chirp + Brp.*Qd.*chir + B.*Qdd.*chir.*chirp + B.*Qd.*chirrp;
Grzp = Br.*Qd.*chizp + B.*Qdd.*chir.*chizp + B.*Qd.*chirzp;
Gzrp = Brp.*Qd.*chiz + B.*Qdd.*chiz.*chirp + B.*Qd.*chizrp;
Gzzp = B.*Qdd.*chiz.*chizp + B.*Qd.*chizzp;
K = (rs/(2*pi)) .* (NRT.*NRS.*Grrp + NRT.*NZS.*Grzp + NZT.*NRS.*Gzrp + NZT.*NZS.*Gzzp);
end

% ===========================================================================
function [F,Fp,Fpp] = F21_and_2derivatives(u,m)
% F = P_{m-1/2}(chi), Fp = dF/dchi, Fpp = d2F/dchi2, with x = (1-chi)/2 = -u/2.
x = -0.5*u;
F   = hyp2f1_series(0.5-m, 0.5+m, 1, x, 80);
Fp  = 0.5*(m^2-0.25) .* hyp2f1_series(1.5-m, 1.5+m, 2, x, 80);
Fpp = -(1/8)*(m^2-0.25)*(9/4-m^2) .* hyp2f1_series(2.5-m, 2.5+m, 3, x, 80);
end

% ===========================================================================
function [Rr,Rp,Rpp] = remainder_and_2derivatives(u,m)
% R, R', R'' in Q_{m-1/2} = -1/2 log(chi-1) F + R.  Analytic limits at u==0, else the
% same Q,Q',Q'' values as the kernel so the split reconstructs it algebraically.
near0 = u == 0;
Rr  = Rm1(m)  * ones(size(u));
Rp  = Rp1(m)  * ones(size(u));
Rpp = Rpp1(m) * ones(size(u));
idx = ~near0;
if any(idx(:))
    ci = u(idx); tv = ci(:).'; maxn = max(m,1);
    [qm,qmd] = chnk.axissymlap2d.qleg_half_miller_vec(tv, maxn);
    q  = reshape(qm(m+1,:),  size(ci));
    qd = reshape(qmd(m+1,:), size(ci));
    chi = 1 + ci;
    qdd = (2*chi.*qd - (m^2-0.25).*q) ./ (-ci.*(2+ci));    % Legendre ODE, 1-chi^2 = -u(u+2)
    [F,Fp,Fpp] = F21_and_2derivatives(ci, m);
    Rr(idx)  = q   + 0.5*log(ci).*F;
    Rp(idx)  = qd  + F./(2*ci) + 0.5*log(ci).*Fp;
    Rpp(idx) = qdd + Fp./ci - F./(2*ci.^2) + 0.5*log(ci).*Fpp;
end
end

% ===========================================================================
function r = Rm1(m)
% R_m(1) = psi(1)-psi(m+1/2)+1/2 log 2, forward recursion (as in the alpha SLP bricks).
R0 = 2.5*log(2); R1 = 2.5*log(2)-2;
if m==0, r = R0; return; end
if m==1, r = R1; return; end
Ra = R0; Rb = R1;
for k = 2:m
    aa = 4*(k-1)/(2*k-1); bb = (2*k-3)/(2*k-1);
    Rc = aa*Rb - bb*Ra; Ra = Rb; Rb = Rc;
end
r = Rb;
end

function rp = Rp1(m)
% R'(1) from the remainder ODE at u=0.
rp = 0.5*(m^2-0.25)*(Rm1(m)+1) + 0.25;
end

function rpp = Rpp1(m)
% R''(1) from the remainder ODE -u(u+2)R''-2(1+u)R'+(m^2-1/4)R = -(u+2)F'-1/2 F
% differentiated once and evaluated at u=0.
Fp1  = 0.5*(m^2-0.25);
Fpp1 = -(1/8)*(m^2-0.25)*(9/4-m^2);
rpp = ( 1.5*Fp1 + 2*Fpp1 + (m^2-0.25-2)*Rp1(m) ) / 4;
end
