function W = LapDLPAxiModalSpecialquad(z, a0, b0, p, gx, gw, Z, Zp, M, iself)
% Per-panel azimuthal-mode close-eval block for the axisymmetric Laplace DLP.
%   W = LapDLPAxiModalSpecialquad(z, a0, b0, p, gx, gw, Z, Zp, M, iself)  -> (M+1) x p
% W(m+1,:) are the panel weights for azimuthal mode m of the DLP (source-normal derivative
% of G).  Split: Cauchy (meridian d_nu, via Dspecialquad) + log (via Sspecialquad) + smooth
% (Gauss).  (Port of dlp_block / kpieces, alpha test_axissymlap_lapdlp_3_bvp.m; cauchy_w
% replaced by Dspecialquad value.)  iself = self-node index, or [] for near-not-self.

hm=(b0-a0)/2; mid=(a0+b0)/2; tt=mid+hm*gx;
y=Z(tt); yp=Zp(tt); spd=abs(yp).'; tau=yp./abs(yp); nv=-1i*tau;     % source tangent/normal (p x 1)
rho=real(z); rpr=real(y).';

[LS,Gm,smo]=kpieces_local(z,y,nv,M);                               % (M+1) x p each
if ~isempty(iself)                                                % closed-form self values
    rp_i=real(y(iself)); nrp_i=real(nv(iself));
    LS(:,iself) =(1/(2*pi))*sqrt(rp_i/rho);
    Gm(:,iself) =-nrp_i/(4*pi*rp_i);
    smo(:,iself)=(nrp_i/(4*pi*rp_i))*(1-Rm1all(M));
end

src.x=y(:); src.wxp=yp(:).*gw(:)*hm; src.nx=nv(:); tgt.x=z;
Aspq=Sspecialquad(tgt,src,Z(a0),Z(b0),'e'); Aspq=Aspq(:).';        % log weights (continuous)
[Acau,~]=Dspecialquad(tgt,src,Z(a0),Z(b0),'e'); Acau=Acau(:).';    % Cauchy value weights
if ~isempty(iself)                                                % self node: principal value
    [Aci,~]=Dspecialquad(tgt,src,Z(a0),Z(b0),'i'); Acau=(Acau+Aci(:).')/2;
end
nvr=nv(:).'; taur=tau(:).'; gwr=gw.';

W=zeros(M+1,p);
for m=0:M
    logblk    = 2*pi*Aspq.*Gm(m+1,:);                             % -log|r-r'| * Gm
    smoothblk = (0.5*log(2*rho.*rpr).*Gm(m+1,:) + smo(m+1,:)).*(spd.*gwr*hm);
    cauchyblk = real( 2*pi*1i*Acau .* (nvr./taur .* LS(m+1,:)) ); % d_nu * LS  (= Re[n_src/(r-r')] LS)
    W(m+1,:)  = cauchyblk + logblk + smoothblk;
end
end

% ===========================================================================
function [LS,Gm,smo]=kpieces_local(z,y,nv,M)
% per-mode DLP split coefficients over panel nodes (port of kpieces, alpha lapdlp).
Mc=max(M,1); rho=real(z); zt=imag(z);
rp=real(y).'; zp=imag(y).'; nrp=real(nv).'; nzp=imag(nv).';        % 1 x p
chim1=((rho-rp).^2+(zt-zp).^2)./(2*rho*rp); rhat2=(rho-rp).^2+(zt-zp).^2;
chim1=max(chim1,1e-15);
[qm,qmd]=chnk.axissymlap2d.qleg_half_miller_vec(chim1,Mc);
n=numel(chim1); qm=reshape(qm,Mc+1,n); qmd=reshape(qmd,Mc+1,n);
F=F21all(chim1,Mc); chi=1+chim1;
Fp=zeros(Mc+1,n); Fp(1,:)=0.5*(F(2,:)-chi.*F(1,:))./(chi.^2-1);
for m=1:Mc, Fp(m+1,:)=(m-0.5)*(chi.*F(m+1,:)-F(m,:))./(chi.^2-1); end
dchi_drp=-(rho-rp)./(rho*rp)-rhat2./(2*rho*rp.^2); dchi_dzp=-(zt-zp)./(rho*rp);
G=nrp.*dchi_drp+nzp.*dchi_dzp;
pref=(1/(4*pi))*2*sqrt(rp/rho); prefp=(1/(4*pi))./sqrt(rp*rho);
KS=pref.*qm;
KD=nrp.*(prefp.*qm+pref.*qmd.*dchi_drp)+nzp.*pref.*qmd.*dchi_dzp - nrp.*KS./rp;
LS=pref.*F;
Gm=(prefp.*nrp-pref.*nrp./rp).*F + pref.*G.*Fp;
dnu=(nrp.*(rho-rp)+nzp.*(zt-zp))./rhat2;
smo=KD - dnu.*LS + 0.5*log(chim1).*Gm;
LS=LS(1:M+1,:); Gm=Gm(1:M+1,:); smo=smo(1:M+1,:);
end

% ===========================================================================
function F = F21all(chim1, M)
F = zeros(M+1, numel(chim1)); F(1,:) = f21_0(chim1);
if M>=1, F(2,:) = f21_1(chim1); end
chi = 1+chim1;
for k=2:M, aa=4*(k-1)/(2*k-1); bb=(2*k-3)/(2*k-1); F(k+1,:) = aa*chi.*F(k,:) - bb*F(k-1,:); end
end
function Rv = Rm1all(M)
Rv = zeros(M+1,1); Rv(1)=2.5*log(2); if M>=1, Rv(2)=2.5*log(2)-2; end
for k=2:M, aa=4*(k-1)/(2*k-1); bb=(2*k-3)/(2*k-1); Rv(k+1) = aa*Rv(k) - bb*Rv(k-1); end
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
