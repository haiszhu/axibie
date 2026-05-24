function W = LapAxiModalSpecialquad(z, a0, b0, p, gx, gw, Z, Zp, M, iself)
% Per-panel azimuthal-mode close-eval block for the axisymmetric Laplace SLP.
%   W = LapAxiModalSpecialquad(z, a0, b0, p, gx, gw, Z, Zp, M, iself)  -> (M+1) x p
% W(m+1,:) are the panel weights for azimuthal mode m, so that
%   (S_m sigma)(z) contribution of this panel = W(m+1,:) * sigma_panel,
% with the physical per-mode SLP kernel  K^S_m = (1/2pi) sqrt(rho'/rho) Q_{m-1/2}(chi).
% Self/near panel: log singularity via Sspecialquad (Helsing/Barnett log product
% integration) + smooth Gauss remainder; far panels use plain Gauss (driver, not here).
%   z      : complex target (rho + i z)
%   a0,b0  : source-panel parameter endpoints
%   p,gx,gw: panel order and Legendre nodes/weights on [-1,1]
%   Z,Zp   : meridian curve and derivative handles
%   M      : max azimuthal mode (returns modes 0..M)
%   iself  : self-node index on this panel, or [] for near-but-not-self
% (Port of nearself_block, flag=true branch, alpha test_axissymlap_lapslp_3_bvp.m.)

hm=(b0-a0)/2; mid=(a0+b0)/2; tt=mid+hm*gx;
y=Z(tt); yp=Zp(tt); rp=real(y); rhoT=real(z);
rhat=abs(z-y); chim1=((rhoT-rp).^2+(imag(z)-imag(y)).^2)./(2*rhoT*rp);
Qall=reshape(chnk.axissymlap2d.qleg_half_miller_vec(chim1(:).',M), M+1, numel(chim1));
Fall=F21all(chim1(:).', M);
Rall=Rm1all(M);
pref=((1/(4*pi))*2.*sqrt(rp./rhoT)).'; spd=abs(yp).'; gwr=gw.'; rhatr=rhat.';

% log product integration (Sspecialquad): A = -1/(2pi) log|z-y| panel weights
src.x=y(:); src.wxp=yp(:).*gw(:)*hm; src.nx=-1i*yp(:)./abs(yp(:)); tgt.x=z;
A=Sspecialquad(tgt,src,Z(a0),Z(b0),'e'); A=A(:).';

W=zeros(M+1,p);
for m=0:M
    F=Fall(m+1,:); Q=Qall(m+1,:);
    smoothQ=Q+log(rhatr).*F;                                  % Q + log|r-r'| F  (smooth)
    if ~isempty(iself), smoothQ(iself)=0.5*log(2*rhoT^2)+Rall(m+1); end
    w_smooth=pref.*smoothQ.*spd.*gwr*hm;
    w_sing  =2*pi*A.*(pref.*F);                               % -log|r-r'| F via Sspecialquad
    W(m+1,:)=w_smooth+w_sing;
end
end

% ===========================================================================
function F = F21all(chim1, M)                       % [M+1 x n]: F21^(m) = P_{m-1/2}(chi)
F = zeros(M+1, numel(chim1));
F(1,:) = f21_0(chim1);
if M>=1, F(2,:) = f21_1(chim1); end
chi = 1+chim1;
for k=2:M
    aa=4*(k-1)/(2*k-1); bb=(2*k-3)/(2*k-1);
    F(k+1,:) = aa*chi.*F(k,:) - bb*F(k-1,:);
end
end
function Rv = Rm1all(M)                             % [M+1 x 1]: R_m(1) self log-remainder
Rv = zeros(M+1,1); Rv(1)=2.5*log(2); if M>=1, Rv(2)=2.5*log(2)-2; end
for k=2:M
    aa=4*(k-1)/(2*k-1); bb=(2*k-3)/(2*k-1);
    Rv(k+1) = aa*Rv(k) - bb*Rv(k-1);
end
end
function f = f21_0(t)
m1 = t./(2+t); f = zeros(size(t)); cn = ones(size(t));
for n=0:80, if n>0, cn = cn.*((2*n-1)/(2*n))^2.*m1; end; f = f + cn; end
f = sqrt(1-m1).*f;
end
function f = f21_1(t)
m1 = t./(2+t); f = zeros(size(t)); cn = ones(size(t));
for n=0:80, if n>0, cn = cn.*((2*n-1)/(2*n))^2.*m1; end; f = f + (1+4*n).*cn; end
f = sqrt(1-m1).*f;
end
