function A = LapDLPGradAxiSpecialMat(t,s)
% DLP prime 

np = s.np; t_in = s.tpan(:);
% refine 1st panel
lam = 1/2; n_split1 = 1;
while t_in(2) > s.t(1)
    split = (1-lam)*t_in(1) + lam*t_in(2);
    t_in = [t_in(1); split; t_in(2:end)];
    n_split1 = n_split1+1;
end
idx1 = 1:s.p*n_split1;
% refine last panel
n_split2 = n_split1+np-2;
while t_in(end-1) < s.t(end)
    split = (1-lam)*t_in(end) + lam*t_in(end-1);
    t_in = [t_in(1:end-1); split; t_in(end)];
end
tpan = t_in;
idx2 = (s.p*n_split2+1):(s.p*(numel(tpan)-1));

% form interpolation matrix...
s_aux = s; s_aux.tpan = tpan; s_aux = quadr(s_aux,[],'p','G');
t_aux1 = 2*(s_aux.t(idx1)-(s.tlo(1)+s.thi(1))/2)/(s.thi(1)-s.tlo(1));
t_aux2 = 2*(s_aux.t(idx2)-(s.tlo(end)+s.thi(end))/2)/(s.thi(end)-s.tlo(end));
L1 = interpmat_1d(t_aux1,gauss(s.p)); L2 = interpmat_1d(t_aux2,gauss(s.p));

% aux system... larger than A
Amat = zeros(numel(t.x),numel(s_aux.x));

for k=1:s_aux.np
    % kth panel...
    pidx = (k-1)*s.p+(1:s.p);
    sk = []; sk.p = s_aux.p; sk.Z = s_aux.Z; sk.tlo = s_aux.tlo(k); sk.thi = s_aux.thi(k);
    sk.np = 1; sk.xlo = s_aux.xlo(k); sk.xhi = s_aux.xhi(k); sk.w = s_aux.w(pidx);
    sk.x = s_aux.x(pidx); sk.xp = s_aux.xp(pidx); sk.xpp = s_aux.xpp(pidx); sk.sp = s_aux.sp(pidx);
    sk.tang = s_aux.tang(pidx); sk.nx = s_aux.nx(pidx); sk.cur = s_aux.cur(pidx); sk.ws = s_aux.ws(pidx);
    sk.t = s_aux.t(pidx); sk.wxp = s_aux.wxp(pidx);

    % close target index...
    panlen = sum(sk.ws); c = 1.5; ik = (abs(t.x - sk.xlo) + abs(t.x - sk.xhi)) < c*panlen;

    % build interaction matrix for close targets...
    if sum(ik)>0
        tk = []; tk.x = t.x(ik(:)); tk.nx = t.nx(ik(:));
        Ak = AxiDGradSpecialquad(tk,sk);
    end

    % naive implementation for far targets interaction... (source normal)
    tfar = []; tfar.x = t.x(~ik(:)); tfar.nx = t.nx(~ik(:));
    x1 = real(tfar.x(:)); x2 = imag(tfar.x(:)); nx1 = real(tfar.nx(:)); nx2 = imag(tfar.nx(:));
    y1 = real(sk.x(:))'; y2 = imag(sk.x(:))'; ny1 = real(sk.nx(:))'; ny2 = imag(sk.nx(:))';
    dxt = bsxfun(@minus,x1,y1); dzt = bsxfun(@minus,x2,y2); r2 = dxt.^2+dzt.^2; rr = bsxfun(@times,x1,y1);
    chi = 1+r2./(2*rr); tmp = 2./(chi+1);
    [SmatK,SmatE] = ellipke(tmp); SnaiveK = sqrt(tmp).*SmatK; SnaiveE = sqrt(tmp).*SmatE;
    sqrty1x1 = sqrt(bsxfun(@rdivide,y1,x1));
    md = bsxfun(@times,nx1,dxt)+bsxfun(@times,nx2,dzt); mpd = bsxfun(@times,ny1,dxt)+bsxfun(@times,ny2,dzt);
    Phi = (md-bsxfun(@rdivide,bsxfun(@times,nx1,r2),2*x1)).*(mpd+bsxfun(@rdivide,bsxfun(@times,ny1,r2),2*y1));
    nn = bsxfun(@times,nx1,ny1); mm = bsxfun(@times,nx2,ny2);
    g = bsxfun(@times,bsxfun(@times,nx2,ny1),x1)-bsxfun(@times,bsxfun(@times,nx1,ny2),y1)-nn.*dzt;
    sK = nn./(2*rr)+Phi./(rr.*(chi+1).*r2);
    sE = 2*(chi.*nn+mm)./r2-8*chi.*Phi./((chi+1).*r2.^2)+3*dzt.*g./(rr.*r2);
    AfarK = 1/(4*pi)*(SnaiveK.*(sqrty1x1.*sK)+SnaiveE.*(sqrty1x1.*sE))*diag(sk.ws);

    % add contribution of kth panel to close and far targets...
    if sum(ik)>0, Amat( ik(:),pidx) = Amat( ik(:),pidx) + Ak; end
    Amat(~ik(:),pidx) = Amat(~ik(:),pidx) + AfarK;

end

intMat = blkdiag(L1,eye((np-2)*s.p),L2);
A = Amat*intMat;

end

% ==============================================================================
function A = AxiDGradSpecialquad(t,s)
% axisymmtric Laplace DLP prime special quadrature on one panel...
% whoever calls this function need to make sure "t.x" is close to panel "s"

be = 2; sf = quadr_panf(s,be,'G');
[sq,coeffs1,coeffs2,coeffs3,coeffs4] = DGradCoeffs(t,sf);
Imat = interpmat(s.p,sf.p,'G');
SspecialMat = Sspecialquad(t,sf,s.xlo,s.xhi,'e');
[Dval,Dz] = Dspecialquad(t,sf,s.xlo,s.xhi,'e');
Alog     = 2*pi*(SspecialMat.*(sq.*coeffs1))*Imat;
Asmooth  = (sq.*coeffs2)*diag(sf.ws)*Imat;
Acauchy  = 2*pi*real(Dval.*bsxfun(@rdivide,sq.*coeffs3,sf.nx.'))*Imat;
Adcauchy = -2*pi*real(Dz.*bsxfun(@rdivide,sq.*coeffs4,sf.nx.'))*Imat;
A = 1/(4*pi)*(Alog + Asmooth + Acauchy + Adcauchy);
end

% ------------------------------------------------------------------------------
function [sq,coeffs1,coeffs2,coeffs3,coeffs4] = DGradCoeffs(t,s)
x1=real(t.x(:)); x2=imag(t.x(:)); nx1=real(t.nx(:)); nx2=imag(t.nx(:));
y1=real(s.x(:)).'; y2=imag(s.x(:)).'; ny1=real(s.nx(:)).'; ny2=imag(s.nx(:)).';
dxt=bsxfun(@minus,x1,y1); dzt=bsxfun(@minus,x2,y2); r2=dxt.^2+dzt.^2; rr=bsxfun(@times,x1,y1);
chi=1+r2./(2*rr); sq=sqrt(bsxfun(@rdivide,y1,x1));
md=bsxfun(@times,nx1,dxt)+bsxfun(@times,nx2,dzt); mpd=bsxfun(@times,ny1,dxt)+bsxfun(@times,ny2,dzt);
Phi=(md-bsxfun(@rdivide,bsxfun(@times,nx1,r2),2*x1)).*(mpd+bsxfun(@rdivide,bsxfun(@times,ny1,r2),2*y1));
nn=bsxfun(@times,nx1,ny1); mm=bsxfun(@times,nx2,ny2);
g=bsxfun(@times,bsxfun(@times,nx2,ny1),x1)-bsxfun(@times,bsxfun(@times,nx1,ny2),y1)-nn.*dzt;
sK=nn./(2*rr)+Phi./(rr.*(chi+1).*r2);
sE=2*(chi.*nn+mm)./r2-8*chi.*Phi./((chi+1).*r2.^2)+3*dzt.*g./(rr.*r2);
F = TruncatedLegendref(t.x,s.x); [TrLe,~] = TruncatedLegendreEf(t.x,s.x);
[KK,EE]=ellipke(2./(chi+1));
Rt=sqrt(2./(chi+1)).*KK+1/2*log(r2).*F;
Rb=sqrt(2./(chi+1)).*EE+1/2*log(r2).*TrLe./(chi+1);
n=nx1+1i*nx2; nu=ny1+1i*ny2;
coeffs1 = sK.*F+sE.*TrLe./(chi+1);
coeffs2 = nn./(2*rr).*Rt+(3*chi.*nn-mm)./(rr.*(chi+1)).*Rb+(Phi-md.*mpd)./(rr.*(chi+1).*r2).*Rt;
coeffs3 = (-4*chi./(chi+1).*(bsxfun(@times,n,ny1./y1)-bsxfun(@times,nx1./x1,nu))+3i*g./rr).*Rb+bsxfun(@times,md,nu).*Rt./(rr.*(chi+1));
coeffs4 = -4*chi./(chi+1).*bsxfun(@times,n,nu).*Rb;
end

% ------------------------------------------------------------------------------
function A = TruncatedLegendref(tx,sx)
x1 = real(tx(:)); x2 = imag(tx(:));
y1 = real(sx(:))'; y2 = imag(sx(:))';
tmp = bsxfun(@rdivide,bsxfun(@minus,x1,y1),x1);
if 0
    T = bsxfun(@rdivide,bsxfun(@minus,x1,y1).^2+bsxfun(@minus,x2,y2).^2,4*x1.^2).*(1+tmp+tmp.^2+tmp.^3+tmp.^4);
    A = 1-1/4*T+9/64*T.^2-25/256*T.^3+1225/16384*T.^4;
else
    T = bsxfun(@rdivide,bsxfun(@minus,x1,y1).^2+bsxfun(@minus,x2,y2).^2,4*x1.^2).*(1+tmp.*(1+tmp.*(1+tmp)));
    A = 1-1/4*T.*(1-9/16*T.*(1-25/36*T));
end
end

% ------------------------------------------------------------------------------
function [A,Al] = TruncatedLegendreEf(tx,sx)
x1 = real(tx(:)); x2 = imag(tx(:));
y1 = real(sx(:))'; y2 = imag(sx(:))';
tmp = bsxfun(@rdivide,bsxfun(@minus,x1,y1),x1);
if 1
    r2 = bsxfun(@minus,x1,y1).^2+bsxfun(@minus,x2,y2).^2;
    T = bsxfun(@rdivide,r2,4*x1.^2).*(1+tmp+tmp.^2+tmp.^3+tmp.^4);
    Al = 1-1/8*T+3/64*T.^2-25/1024*T.^3;
    A = T.*Al;
    Al = A./r2;
else
    T = bsxfun(@rdivide,bsxfun(@minus,x1,y1).^2+bsxfun(@minus,x2,y2).^2,4*x1.^2).*(1+tmp.*(1+tmp.*(1+tmp)));
    A = 1-1/4*T.*(1-9/16*T.*(1-25/36*T));
end
end
