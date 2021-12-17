function [A,s_aux]= AxiStokesSpecialMat(t,s)
% special quadrature matrix....
%
% also return s_aux for test purpose...
%
%
% Hai 11/27/20

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
%     n_split2 = n_split2+1;
end
tpan = t_in;
idx2 = (s.p*n_split2+1):(s.p*(numel(tpan)-1));

% form interpolation matrix...
s_aux = s; s_aux.tpan = tpan; s_aux = quadr(s_aux,[],'p','G'); % maybe could avoid regenerating via interpolation...
t_aux1 = 2*(s_aux.t(idx1)-(s.tlo(1)+s.thi(1))/2)/(s.thi(1)-s.tlo(1)); 
t_aux2 = 2*(s_aux.t(idx2)-(s.tlo(end)+s.thi(end))/2)/(s.thi(end)-s.tlo(end)); 
L1 = interpmat_1d(t_aux1,gauss(s.p)); L2 = interpmat_1d(t_aux2,gauss(s.p));

% aux system... larger than A
Amat = zeros(numel(t.x),2*numel(s_aux.x));

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
    tk = []; tk.x = t.x(ik(:)); 
    
    % build interaction matrix for close targets...
    Ak = AxiStokesSpecialquad(tk,sk);    
    
    % naive implementation for far targets interaction...
    tfar = []; tfar.x = t.x(~ik(:));
    x1 = real(tfar.x(:)); x2 = imag(tfar.x(:));     % target (colume)
    y1 = real(sk.x(:))'; y2 = imag(sk.x(:))';   % source (row)
    r2 = bsxfun(@minus,x1,y1).^2+bsxfun(@minus,x2,y2).^2;
    chi = 1+r2./(2*bsxfun(@times,x1,y1)); tmp = 2./(chi+1);
    [SmatK,SmatE] = ellipke(tmp);
    SnaiveK = sqrt(tmp).*SmatK; SnaiveE = sqrt(tmp).*SmatE;
    sqrty1x1 = sqrt(bsxfun(@rdivide,y1,x1));
    S11 = sqrty1x1.*(bsxfun(@plus,x1.^2,y1.^2)./bsxfun(@times,x1,y1)+2*bsxfun(@minus,y2,x2).^2./bsxfun(@times,x1,y1));
    S12 = sqrty1x1.*bsxfun(@rdivide,bsxfun(@minus,x2,y2),x1);
    S21 = sqrty1x1.*bsxfun(@rdivide,bsxfun(@minus,y2,x2),y1);
    S22 = sqrty1x1*2;
    AfarK = 1/(8*pi)*[(SnaiveK.*S11)*diag(sk.ws),(SnaiveK.*S12)*diag(sk.ws);(SnaiveK.*S21)*diag(sk.ws),(SnaiveK.*S22)*diag(sk.ws)];
    S11E = sqrty1x1.*(-2-4*chi+2*chi.*bsxfun(@minus,y1,x1).^2./r2);
    S12E = 2*sqrty1x1.*(bsxfun(@minus,x1,y1).*bsxfun(@minus,x2,y2)./r2-bsxfun(@rdivide,bsxfun(@minus,x2,y2),2*x1));
    S21E = 2*sqrty1x1.*(bsxfun(@minus,x1,y1).*bsxfun(@minus,x2,y2)./r2+bsxfun(@rdivide,bsxfun(@minus,x2,y2),2*y1));
    S22E = 2*sqrty1x1.*bsxfun(@minus,x2,y2).^2./r2;
    AfarE = 1/(8*pi)*[(SnaiveE.*S11E)*diag(sk.ws),(SnaiveE.*S12E)*diag(sk.ws);(SnaiveE.*S21E)*diag(sk.ws),(SnaiveE.*S22E)*diag(sk.ws)];
    
    % add contribution of kth panel to close and far targets...
    Amat(ik(:),[pidx,end/2+pidx]) = Amat(ik(:),[pidx,end/2+pidx]) + Ak(1:end/2,:) + 1i*Ak(end/2+1:end,:);
    Amat(~ik(:),[pidx,end/2+pidx]) = Amat(~ik(:),[pidx,end/2+pidx]) + AfarK(1:end/2,:) + 1i*AfarK(end/2+1:end,:) + AfarE(1:end/2,:) + 1i*AfarE(end/2+1:end,:);
    
    
end

Aaux = [real(Amat);imag(Amat)];
intMat = blkdiag(L1,eye((np-2)*s.p),L2);
A = Aaux*blkdiag(intMat,intMat);

% keyboard

end

function A = AxiStokesSpecialquad(t,s)
% axisymmtric Stokes special quadrature on one panel...
%
% whoever calls this function need to make sure "t.x" is close to panel "s"
%
% Hai 11/24/20

Asingular = AxiStokesSingularquad(t,s);
Asmooth = AxiStokesSmoothquad(t,s);
% A = 1/(8*pi)*(Asingular + Asmooth);
AsmoothE = AxiStokesSmoothEquad(t,s);
A = 1/(8*pi)*(Asingular + Asmooth + AsmoothE);

end

function A = AxiStokesSingularquad(t,s)
% special quadrature for close evaluation of axisymmetric stokes kernel... 
% on one panel...
%
% t.x: close target, s.x: source; also need s.xlo, s.xhi, etc... 
%
% Hai 11/23/20

be = 2; % upsampling factor...
% singular part...
sf = quadr_panf(s,be,'G'); % upsampled panel
x1 = real(t.x(:)); x2 = imag(t.x(:));     % source (colume)
y1 = real(sf.x(:))'; y2 = imag(sf.x(:))';   % target (row)
TrLe = TruncatedLegendref(x1+1i*x2,y1+1i*y2);
sqrty1x1 = sqrt(bsxfun(@rdivide,y1,x1));
S11 = sqrty1x1.*(bsxfun(@plus,x1.^2,y1.^2)./bsxfun(@times,x1,y1)+2*bsxfun(@minus,y2,x2).^2./bsxfun(@times,x1,y1));
S12 = sqrty1x1.*bsxfun(@rdivide,bsxfun(@minus,x2,y2),x1);
S21 = sqrty1x1.*bsxfun(@rdivide,bsxfun(@minus,y2,x2),y1);
S22 = sqrty1x1*2;
SspecialMat = Sspecialquad(t,sf,s.xlo,s.xhi,'e');
A11 = 2*pi*(SspecialMat.*(TrLe.*S11))*interpmat(s.p,sf.p,'G');
A12 = 2*pi*(SspecialMat.*(TrLe.*S12))*interpmat(s.p,sf.p,'G');
A21 = 2*pi*(SspecialMat.*(TrLe.*S21))*interpmat(s.p,sf.p,'G');
A22 = 2*pi*(SspecialMat.*(TrLe.*S22))*interpmat(s.p,sf.p,'G');
A = [A11,A12;A21,A22];

end

function A = AxiStokesSmoothquad(t,s)
% quadrature for the smooth part of axisymmetric stokes kernel associated with EllipticK part... 
% on one panel...
%
% t.x: close target, s.x: source; also need s.xlo, s.xhi, etc... 
%
% Hai 11/23/20

be = 2; % upsampling factor...
% smooth part...
sf = quadr_panf(s,be,'G'); 
x1 = real(t.x(:)); x2 = imag(t.x(:));     % source (colume)
y1 = real(sf.x(:))'; y2 = imag(sf.x(:))';   % target (row)
r2 = bsxfun(@minus,x1,y1).^2+bsxfun(@minus,x2,y2).^2;
chi = 1+r2./(2*bsxfun(@times,x1,y1)); tmp = 2./(chi+1);
Smat1 = sqrt(tmp).*ellipke(tmp);
Smat2 = 1/2*log(r2).*TruncatedLegendref(x1+1i*x2,y1+1i*y2);
sqrty1x1 = sqrt(bsxfun(@rdivide,y1,x1));
S11 = sqrty1x1.*(bsxfun(@plus,x1.^2,y1.^2)./bsxfun(@times,x1,y1)+2*bsxfun(@minus,y2,x2).^2./bsxfun(@times,x1,y1));
S12 = sqrty1x1.*bsxfun(@rdivide,bsxfun(@minus,x2,y2),x1);
S21 = sqrty1x1.*bsxfun(@rdivide,bsxfun(@minus,y2,x2),y1);
S22 = sqrty1x1*2;
A11 = ((Smat1+Smat2).*S11)*diag(sf.ws)*interpmat(s.p,sf.p,'G');
A12 = ((Smat1+Smat2).*S12)*diag(sf.ws)*interpmat(s.p,sf.p,'G');
A21 = ((Smat1+Smat2).*S21)*diag(sf.ws)*interpmat(s.p,sf.p,'G');
A22 = ((Smat1+Smat2).*S22)*diag(sf.ws)*interpmat(s.p,sf.p,'G');
A = [A11,A12;A21,A22];

end

function A = AxiStokesSmoothEquad(t,s)
% quadrature for the smooth part of axisymmetric stokes kernel associated with EllipticE part... 
% on one panel...
% nothing special...
% t.x: close target, s.x: source; also need s.xlo, s.xhi, etc... 
%
% Hai 11/23/20

be = 2; % upsampling factor...
% smooth part...
sf = quadr_panf(s,be,'G'); 
x1 = real(t.x(:)); x2 = imag(t.x(:));     % source (colume)
y1 = real(sf.x(:))'; y2 = imag(sf.x(:))';   % target (row)
r2 = bsxfun(@minus,x1,y1).^2+bsxfun(@minus,x2,y2).^2;
chi = 1+r2./(2*bsxfun(@times,x1,y1)); tmp = 2./(chi+1);
[~,Smat1] = ellipke(tmp); Smat1 = sqrt(tmp).*Smat1; % take the elliptic integral of second kind...
sqrty1x1 = sqrt(bsxfun(@rdivide,y1,x1)); x1y1 = bsxfun(@times,x1,y1);
SspecialMat = Sspecialquad(t,sf,s.xlo,s.xhi,'e');
DspecialMat = Dspecialquad(t,sf,s.xlo,s.xhi,'e');
[TrLe,TrLel] = TruncatedLegendreEf(x1+1i*x2,y1+1i*y2);
Smat2 = 1/2*log(r2).*TrLe;
Smat11 = sqrty1x1.*((-2-4*chi).*r2+2*chi.*bsxfun(@minus,x1,y1).^2)./(chi+1).*TrLel;
Smat12 = 2*sqrty1x1.*(bsxfun(@minus,x1,y1).*bsxfun(@minus,x2,y2)-bsxfun(@rdivide,bsxfun(@minus,x2,y2),2*x1).*r2)./(chi+1).*TrLel;
Smat21 = 2*sqrty1x1.*(bsxfun(@minus,x1,y1).*bsxfun(@minus,x2,y2)+bsxfun(@rdivide,bsxfun(@minus,x2,y2),2*y1).*r2)./(chi+1).*TrLel;
Smat22 = 2*sqrty1x1.*bsxfun(@minus,x2,y2).^2./(chi+1).*TrLel;
Rhat = Smat1+Smat2./(chi+1);

% A11
A11smooth = sqrty1x1.*(-2-4*chi).*Rhat*diag(sf.ws)*interpmat(s.p,sf.p,'G');
A11cauchy = 4*pi*(DspecialMat.*bsxfun(@rdivide,bsxfun(@minus,x1,y1),sf.nx.')).*(sqrty1x1.*chi).*Rhat;
A11cauchy = real(A11cauchy)*interpmat(s.p,sf.p,'G');
A11log = 2*pi*(SspecialMat.*Smat11)*interpmat(s.p,sf.p,'G');

% A12
A12smooth = -sqrty1x1.*bsxfun(@rdivide,bsxfun(@minus,x2,y2),x1).*Rhat*diag(sf.ws)*interpmat(s.p,sf.p,'G');
A12cauchy = 4*pi*(DspecialMat.*bsxfun(@rdivide,bsxfun(@minus,x2,y2),sf.nx.')).*sqrty1x1.*Rhat;
A12cauchy = real(A12cauchy)*interpmat(s.p,sf.p,'G');
A12log = 2*pi*(SspecialMat.*Smat12)*interpmat(s.p,sf.p,'G');

% A21
A21smooth = sqrty1x1.*bsxfun(@rdivide,bsxfun(@minus,x2,y2),y1).*Rhat*diag(sf.ws)*interpmat(s.p,sf.p,'G');
A21cauchy = A12cauchy;
A21log = 2*pi*(SspecialMat.*Smat21)*interpmat(s.p,sf.p,'G');

% A22
A22smooth = zeros(size(A11smooth));
A22cauchy = 4*pi*(DspecialMat.*bsxfun(@rdivide,1i*bsxfun(@minus,x2,y2),sf.nx.')).*sqrty1x1.*Rhat;
A22cauchy = real(A22cauchy)*interpmat(s.p,sf.p,'G');
A22log = 2*pi*(SspecialMat.*Smat22)*interpmat(s.p,sf.p,'G');

A = [A11smooth+A11cauchy+A11log, A12smooth+A12cauchy+A12log; ...
     A21smooth+A21cauchy+A21log, A22smooth+A22cauchy+A22log];


end

function A = TruncatedLegendref(tx,sx)
% truncated Legendre function...
%
% tx: target; sx: source
%
% Hai 11/20/20

% taylor expansion
% syms x
% taylor(hypergeom([1/2 1/2],1,x), x)

x1 = real(tx(:)); x2 = imag(tx(:));     % source (colume)
y1 = real(sx(:))'; y2 = imag(sx(:))';   % target (row)
tmp = bsxfun(@rdivide,bsxfun(@minus,x1,y1),x1);
if 0 % a 4th order expansion...
    T = bsxfun(@rdivide,bsxfun(@minus,x1,y1).^2+bsxfun(@minus,x2,y2).^2,4*x1.^2).*(1+tmp+tmp.^2+tmp.^3+tmp.^4);
%     A = 1-1/4*T+1/2^6*T.^2-1/(3*2^5)*T.^3+25/2^14*T.^4;
    A = 1-1/4*T+9/64*T.^2-25/256*T.^3+1225/16384*T.^4;
else % a 3rd order expansion... (this one is probably more balanced...)
    T = bsxfun(@rdivide,bsxfun(@minus,x1,y1).^2+bsxfun(@minus,x2,y2).^2,4*x1.^2).*(1+tmp.*(1+tmp.*(1+tmp)));
    A = 1-1/4*T.*(1-9/16*T.*(1-25/36*T));
end


end

function [A,Al] = TruncatedLegendreEf(tx,sx)
% truncated Legendre function for EllipticE...
%
% tx: target; sx: source
%
% A for the smooth part, and Al is for the log associated part...
%
% Hai 11/26/20

% taylor expansion
% syms x
% taylor(hypergeom([-1/2 3/2],1,x), x)

x1 = real(tx(:)); x2 = imag(tx(:));     % source (colume)
y1 = real(sx(:))'; y2 = imag(sx(:))';   % target (row)
tmp = bsxfun(@rdivide,bsxfun(@minus,x1,y1),x1);
if 1 % a 4th order expansion...
    r2 = bsxfun(@minus,x1,y1).^2+bsxfun(@minus,x2,y2).^2;
    T = bsxfun(@rdivide,r2,4*x1.^2).*(1+tmp+tmp.^2+tmp.^3+tmp.^4);
%     A = 1-1/4*T+1/2^6*T.^2-1/(3*2^5)*T.^3+25/2^14*T.^4;
    Al = 1-1/8*T+3/64*T.^2-25/1024*T.^3;
    A = T.*Al;
    Al = A./r2;
else % a 3rd order expansion... (this one is probably more balanced...)
    T = bsxfun(@rdivide,bsxfun(@minus,x1,y1).^2+bsxfun(@minus,x2,y2).^2,4*x1.^2).*(1+tmp.*(1+tmp.*(1+tmp)));
    A = 1-1/4*T.*(1-9/16*T.*(1-25/36*T));
end


end