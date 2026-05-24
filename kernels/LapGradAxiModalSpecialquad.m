function W = LapGradAxiModalSpecialquad(z, ntn, a0, b0, p, gx, gw, Z, Zp, M, iself)
% Per-panel azimuthal-mode close-eval block for the axisymmetric Laplace SLPn S'
% (target-normal derivative of the SLP).
%   W = LapGradAxiModalSpecialquad(z, ntn, a0, b0, p, gx, gw, Z, Zp, M, iself)  -> (M+1) x p
% ntn = target meridian unit normal [n_rho; n_z].  Split: log (Sspecialquad) + target-normal
% Cauchy (Dspecialquad value) + smooth (Gauss).  Per-mode split coefficients (G1,Gc,G0)
% come from the validated util axissymlap_slpn_kernel_split (chunkie convention), rescaled
% to the physical convention (1/(-2pi rho_src)) used by the other gamma routines.
% iself = self-node index, or [] for near-but-not-self.

hm=(b0-a0)/2; mid=(a0+b0)/2; tt=mid+hm*gx;
y=Z(tt); yp=Zp(tt); spd=abs(yp).'; tau=yp./abs(yp); nv=-1i*tau;     % source tangent/normal
rp=real(y).'; rho=real(z); nt_c=ntn(1)+1i*ntn(2);
src=struct('r',[real(y(:))';imag(y(:))'],'n',[real(nv(:))';imag(nv(:))']);
tgt=struct('r',[real(z);imag(z)],'n',ntn(:));

ss.x=y(:); ss.wxp=yp(:).*gw(:)*hm; ss.nx=nv(:); tg.x=z;
Aspq=Sspecialquad(tg,ss,Z(a0),Z(b0),'e'); Aspq=Aspq(:).';
[Acau,~]=Dspecialquad(tg,ss,Z(a0),Z(b0),'e'); Acau=Acau(:).';
if ~isempty(iself)                                                 % self node: principal value
    [Aci,~]=Dspecialquad(tg,ss,Z(a0),Z(b0),'i'); Acau=(Acau+Aci(:).')/2;
end
taur=tau(:).'; gwr=gw.'; convp=1./(-2*pi*rp);                       % chunkie -> physical, per node

W=zeros(M+1,p);
for m=0:M
    [G1,Gc,G0]=axissymlap_slpn_kernel_split(src, tgt, m);          % chunkie conv, 1 x p
    G1p=G1.*convp; Gcp=Gc.*convp; G0p=G0.*convp;
    logblk    = -2*pi*Aspq.*G1p;
    cauchyblk = real(2*pi*1i*nt_c*Acau.*(Gcp./taur));             % ((r-r').n_t)/rhat^2 * Gc
    smoothblk = G0p.*(spd.*gwr*hm);
    W(m+1,:)  = logblk + cauchyblk + smoothblk;
end
end

function [G1, Gc, G0] = axissymlap_slpn_kernel_split(srcinfo, targinfo, m)
% Split the per-mode axisymmetric Laplace SLPn kernel into
%
%   K = log(|r-r'|) G1 + ((r-r') . n_t)/|r-r'|^2 Gc + G0 .
%
% This matches axissymlap_slpn_kernel's chunkie convention: the returned kernel
% already includes the azimuthal factor and source-radius Jacobian.  Here n_t is
% the target meridian normal.  G1, Gc, and G0 are bounded at coincidence; the
% explicit Cauchy factor carries the 1/|r-r'| singularity.

src_r  = srcinfo.r(1,:);   src_z  = srcinfo.r(2,:);
targ_r = targinfo.r(1,:);  targ_z = targinfo.r(2,:);
nr = targinfo.n(1,:).';    nz = targinfo.n(2,:).';

ns = numel(src_r);  nt = numel(targ_r);
rt = repmat(targ_r.',1,ns);  zt = repmat(targ_z.',1,ns);
rs = repmat(src_r,   nt,1);  zs = repmat(src_z,   nt,1);
nr = repmat(nr,1,ns);        nz = repmat(nz,1,ns);

dr = rt - rs;
dz = zt - zs;
rhat2 = dr.^2 + dz.^2;
chim1 = rhat2 ./ (2*rt.*rs);
chi = 1 + chim1;

pref = (1/(2*pi)) .* sqrt(rs./rt) ./ (rs.*rt);
A = pref .* nr .* rs/2;
B = nr .* (rt - chi.*rs) + nz .* dz;
C = -pref .* B;

[F, Fp] = F21_and_derivative(chim1, m);
L = A.*F + C.*Fp;
G1 = -L;
Gc = pref .* F .* rt .* rs;

[R, Rp] = remainder_and_derivative(chim1, m);

% The q' singular term contains B/(2*(chi-1)).  Pulling out the explicit
% Cauchy factor leaves this bounded correction in the smooth remainder.
smooth_from_cauchy = -pref .* F .* nr .* rs/2;
G0 = 0.5*log(2*rt.*rs).*L + smooth_from_cauchy + A.*R + C.*Rp;

end

% ===========================================================================
function [F,Fp] = F21_and_derivative(chim1,m)
% F = P_{m-1/2}(chi) and Fp = dF/dchi from the hypergeometric coefficient in
% the Helsing split, with x=(1-chi)/2=-chim1/2.
x = -0.5*chim1;
F = hyp2f1_series(0.5-m, 0.5+m, 1, x, 80);
Fp = 0.5*(m^2-0.25) .* hyp2f1_series(1.5-m, 1.5+m, 2, x, 80);
end

% ===========================================================================
function [R,Rp] = remainder_and_derivative(chim1,m)
% R and dR/dchi in Q_{m-1/2} = -1/2 log(chi-1) F + R.  At exact self nodes use
% analytic limits; for nonzero separations keep the same Q/Q' values as the full
% kernel so the split reconstructs it algebraically.
near0 = chim1 == 0;
R = Rm1(m) * ones(size(chim1));
Rp = Rp1(m) * ones(size(chim1));

idx = ~near0;
if any(idx(:))
    ci = chim1(idx);                                  % preserves chim1 orientation
    tv = ci(:).';
    maxn = max(m,1);
    [qm,qmd] = chnk.axissymlap2d.qleg_half_miller_vec(tv, maxn);
    q  = reshape(qm(m+1,:),  size(ci));
    qd = reshape(qmd(m+1,:), size(ci));
    [F,Fp] = F21_and_derivative(ci, m);
    R(idx)  = q  + 0.5*log(ci).*F;
    Rp(idx) = qd + F./(2*ci) + 0.5*log(ci).*Fp;
end
end

% ===========================================================================
function r = Rm1(m)
% R_m(1) = psi(1)-psi(m+1/2)+1/2 log 2, evaluated by the same forward
% recursion used in the alpha Laplace SLP bricks.
R0 = 2.5*log(2);
R1 = 2.5*log(2)-2;
if m==0, r = R0; return; end
if m==1, r = R1; return; end
Ra = R0; Rb = R1;
for k = 2:m
    aa = 4*(k-1)/(2*k-1);
    bb = (2*k-3)/(2*k-1);
    Rc = aa*Rb - bb*Ra;
    Ra = Rb; Rb = Rc;
end
r = Rb;
end

function rp = Rp1(m)
% R'_m(1), from the Legendre remainder ODE after inserting the log split.
r = Rm1(m);
rp = 0.5*(m^2-0.25)*(r+1) + 0.25;
end

% ===========================================================================
function s = hyp2f1_series(a,b,c,x,N)
s = ones(size(x));
term = ones(size(x));
for n = 0:N-1
    term = term .* ((a+n).*(b+n)./((c+n)*(n+1))) .* x;
    s = s + term;
end
end
