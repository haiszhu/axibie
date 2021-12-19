function A = LapDLPSpecialMatDG(t,s,flag)
% special quadrature for Laplace BVP...
% close-far criteria could be modified... 
% 
% need to use on-surface flag, and compute an average... 
%
% Hai 05/21/21

if nargin < 3, flag = 0; end

be = 2; qntype='G'; 
A = LapDLPmat(t,s);
ls=0; len=s.p; np = s.np;
panlen = zeros(1,np); intMat = interpmat(s.p,be*s.p,qntype);
[~, wf] = gauss(be*s.p); 
for k=1:np
    j = (k-1)*s.p+(1:s.p); 
    panlen(k) = sum(s.ws(j)); 
    ik = (abs(t.x - s.xlo(k)) + abs(t.x - s.xhi(k))) < 1.6*panlen(k);
    if sum(ik)>0
        
        wxpf = wf.*((s.thi(k)-s.tlo(k))/2*intMat*s.xp(j));

        zsc = (s.xhi(k)-s.xlo(k))/2; zmid = (s.xhi(k)+s.xlo(k))/2; % rescaling factor and midpoint of src segment
        tt.x=(t.x(ik(:))-zmid)/zsc;
        ss = {}; ss.Z = @(t) (s.Z((s.tlo(k)*(1-t)+s.thi(k)*(1+t))/2)-zmid)/zsc;
        ss.p = be*s.p; ss.tlo = -1; ss.thi = 1; ss = quadr(ss,ss.p,'p','G');
        
        if flag
            % need to separate target into self and close, for self, need to do an average, for close, do as usual
            ikc = ik; ikc(j) = 0; iks = ik&(~ikc);
            ts.x=(t.x(iks(:))-zmid)/zsc; tc.x=(t.x(ikc(:))-zmid)/zsc;
            Flag.flag = 1; Flag.sx = (s.x(j)-zmid)/zsc; % need to know which branch cut, it seems like more related to positive or negative of imag(s.x)
            Aks = D0formChebyquad(ts,ss,-1,1,Flag);
            Akc = D0formChebyquad(tc,ss,-1,1);
            A(iks,ls+1:ls+len) = real(Aks)*intMat;
            A(ikc,ls+1:ls+len) = real(Akc)*intMat;
        else
            % 
            Ak = D0formChebyquad(tt,ss,-1,1);
            A(ik,ls+1:ls+len) = real(Ak)*intMat;
        end
        
    end
    ls=ls+len;
end
end


function [A,Az,Azz] = D0formChebyquad(t,s,a,b,Flag,interpMat,w0f)
% special quadrature based on a differential 0-form approach 
% with cheby basis
%
% Hai 04/26/21

N = numel(t.x);                          % # of targets
p = numel(s.x);                          % assume panel order is # nodes
if nargin<5, flag = 0; else, flag = Flag.flag; end
if nargin<7
    [t0,~] = gauss(p); be = 2; [t0f,w0f] = gauss(be*p);
    tcen = [-3/4,0,3/4]; tlen = [1/4,1/2,1/4]; t0f=t0f*tlen+ones(be*p,1)*tcen; w0f=w0f*tlen;
    t0f=t0f(:); w0f=w0f(:); interpMat = interpmat_1d(t0f,t0);
end
zsc = (b-a)/2; zmid = (b+a)/2; % rescaling factor and midpoint of src segment
y = (s.x-zmid)/zsc; x = (t.x-zmid)/zsc; yxp = s.xp/zsc; % transformed src nodes, targ pts
VT = ones(p,p); VT(:,2) = y(:); for k=3:p, VT(:,k) = 2*y(:).*VT(:,k-1)-VT(:,k-2); end  % Chebyshev mat @ nodes

idxfr = (abs(x+1) + abs(x-1)) > 1.3*2;   % far
idxnr = (abs(x+1) + abs(x-1)) <= 1.3*2;  % near
xfr = x(idxfr); xnr = x(idxnr);

% separate near targets into 2 catogeries
idx1 = (abs(xnr+1) + abs(xnr-1)) <= 1.3*2;  % closest
idx2 = ~idx1;    % middle

% closest 
x1 = xnr(idx1); N1 = numel(x1);
T1 = zeros(p,N1); U1 = zeros(p,N1);
if flag
    p1i = -atan2(-imag(x1),(real(x1).^2+imag(x1).^2-real(x1)))+atan2(imag(x1),(real(x1).^2+imag(x1).^2+real(x1)));
    p1r = 1/2*log(((1-real(x1)).^2+imag(x1).^2)./((1+real(x1)).^2+imag(x1).^2));
    sxNegIdx = imag(Flag.sx)<0; sxPosIdx = imag(Flag.sx)>0;   % for the simplest case, assume only changes once
    p1i(sxNegIdx) = p1i(sxNegIdx) + pi; p1i(sxPosIdx) = p1i(sxPosIdx) - pi; 
%     if sum(s.cur)<0
%         p1i = p1i - pi;
%     else
%         p1i = p1i + pi;
%     end
else
    [IN,ON] = inpolygon(real(x1),imag(x1),[-1;real(y);1],[0;imag(y);0]); idx = IN|ON; % if within another branch
    p1i = -atan2(-imag(x1),(real(x1).^2+imag(x1).^2-real(x1)))+atan2(imag(x1),(real(x1).^2+imag(x1).^2+real(x1)));
    p1r = 1/2*log(((1-real(x1)).^2+imag(x1).^2)./((1+real(x1)).^2+imag(x1).^2));
    if sum(s.cur)<0
        bidx0 = abs(x1)<1e-15;
        p1i(idx) = -2*pi+p1i(idx); p1i(bidx0) = -pi;
    else
        bidx2 = real(x1)<1 & real(x1)>-1 & abs(imag(x1))<1e-15; bidx0 = abs(x1)<1e-15;
        p1i(idx&~bidx2) = 2*pi+p1i(idx&~bidx2); p1i(bidx0) = pi;
    end
end
T1(1,:) = p1r'+1i*p1i'; T1(2,:) = x1.'.*T1(1,:) + 2; k=2; T1(k+1,:) = 2*x1.'.*T1(k,:) - T1(k-1,:);
U1(1,:) = T1(1,:); U1(2,:) = 2*T1(2,:); k=2; U1(k+1,:) = 2*x1.'.*U1(k,:) - U1(k-1,:);
for k=3:p-1
    T1(k+1,:) = (1/k-1/(k-2))-((-1)^k/k-(-1)^(k-2)/(k-2)) + 2*x1.'.*T1(k,:) - T1(k-1,:);
    U1(k+1,:) = 2/k*(1-(-1)^k) + 2*x1.'.*U1(k,:) - U1(k-1,:);
end

% middle (forward recursion fails, instead a solving process)
x2 = xnr(idx2); N2 = numel(x2);
T2 = zeros(p,N2); U2 = zeros(p,N2);
p2i = -atan2(-imag(x2),(real(x2).^2+imag(x2).^2-real(x2)))+atan2(imag(x2),(real(x2).^2+imag(x2).^2+real(x2)));
p2r = 1/2*log(((1-real(x2)).^2+imag(x2).^2)./((1+real(x2)).^2+imag(x2).^2));
T2(1,:) = p2r'+1i*p2i'; T2(2,:) = x2.'.*T2(1,:) + 2; 
U2(1,:) = T2(1,:); U2(2,:) = 2*T2(2,:);
yf = interpMat*y; wxp = w0f.*(interpMat*yxp); 
tmp1 = ones(size(yf)); tmp2 = yf; yfp = zeros(size(yf));
for k=3:p, yfp = 2*yf(:).*tmp2-tmp1; tmp1 = tmp2; tmp2 = yfp; end  % Chebyshev 
tmp1 = ones(size(yf)); tmp2 = 2*yf; yfpU = zeros(size(yf));
for k=3:p, yfpU = 2*yf(:).*tmp2-tmp1; tmp1 = tmp2; tmp2 = yfpU; end  % Chebyshev 2nd kind 
Tp = (-1./bsxfun(@minus,x2,yf.')*(wxp.*yfp)).'; Up = (-1./bsxfun(@minus,x2,yf.')*(wxp.*yfpU)).';
rhs0 = -(1./(3:p-1)-1./(1:p-3)) + ((-1).^(3:p-1)./(3:p-1)-(-1).^(1:p-3)./(1:p-3)); rhs0 = rhs0(:); 
Rhs = rhs0*ones(1,N2); Rhs(1,:) = Rhs(1,:) + T2(2,:); Rhs(end,:) = Rhs(end,:) + Tp; % Rhs = Rhs(:);
h=1/(p-2); E = sin((1:(p-3))'*(1:(p-3))*pi*h); Ei = E\Rhs; Lam = bsxfun(@plus,-2*cos((1:(p-3))'*pi*h),2*x2.'); 
T2(3:p-1,:) = E*(Ei./Lam); T2(end,:) = Tp;
rhsU0 = -2./(3:p-1) + 2*(-1).^(3:p-1)./(3:p-1); rhsU0 = rhsU0(:);
RhsU = rhsU0*ones(1,N2); RhsU(1,:) = RhsU(1,:) + U2(2,:); RhsU(end,:) = RhsU(end,:) + Up; EiU = E\RhsU;
U2(3:p-1,:) = E*(EiU./Lam); U2(end,:) = Up;

% for j=1:numel(x2)
%     Lam = 2*x2(j)-2*cos((1:(p-3))'*pi*h);
%     tmp = E*((1./Lam).*Ei(:,j));
%     T2(3:p-1,j) = tmp.';
% end

% matrix
T = zeros(p,N1+N2); T(:,idx1) = T1; T(:,idx2) = T2; 
Anr = (VT.'\T).'*(1i/(2*pi)); rfr = bsxfun(@minus,xfr,yf.'); wxpfr = bsxfun(@times,wxp,interpMat);
% Afr = 1i/(2*pi)*bsxfun(@times,-1./bsxfun(@minus,xfr,yf.'),wxp.')*interpMat;
Afr = -1i/(2*pi)./rfr*wxpfr; % upsampling for far
A = zeros(N,p); A(idxnr,:) = Anr; A(idxfr,:) = Afr;

if nargout > 1
    % derivative matrix
    U = zeros(p,N1+N2); U(:,idx1) = U1; U(:,idx2) = U2; 
    R =  -(kron(ones(p,1),1./(1-xnr.')) + kron((-1).^(0:p-1).',1./(1+xnr.'))) +...
            repmat((0:p-1)',[1 N1+N2]).*[zeros(1,N1+N2); U(1:p-1,:)];  % hypersingular kernel weights
    Aznr = (VT.'\R).'*(1i/(2*pi*zsc));  % solve for targ complex-deriv mat & rescale
    Azfr = 1i/(2*pi)./rfr.^2*wxpfr; % upsampling for far
    Az = zeros(N,p); Az(idxnr,:) = Aznr; Az(idxfr,:) = Azfr;

    % 2nd derivative matrix...
    RU = zeros(size(R)); RU(1,:) = R(1,:); RU(2,:) = 2*R(2,:);
    for k=3:p, RU(k,:) = RU(k-2,:)+2*R(k,:); end
    S = -(kron(ones(p,1),1./(1-xnr.').^2) - kron((-1).^(0:p-1).',1./(1+xnr.').^2))/2 +...
           repmat((0:p-1)',[1 N1+N2]).*[zeros(1,N1+N2); RU(1:p-1,:)]/2; % supersingular kernel weights
    Azznr = (VT.'\S).'*(1i/(2*pi*zsc^2));
    Azzfr = -1i/(2*pi)./rfr.^3*wxpfr; % upsampling for far
    Azz = zeros(N,p); Azz(idxnr,:) = Azznr; Azz(idxfr,:) = Azzfr;
end

% keyboard

end

function P = interpmat(n,m, qntype) % interpolation matrix from n-pt to m-pt Gauss nodes
% INTERPMAT - create interpolation matrix from n-pt to m-pt Gauss nodes
%
% P = interpmat(n,m) returns a m*n matrix which maps func values on n-pt Gauss-
% Legendre nodes on [-1,1] to values on m-pt nodes.
% Does it the Helsing way via backwards-stable ill-cond Vandermonde solve.
if m==n, P = eye(n); return, end
if qntype=='G', x = gauss(n); y = gauss(m); 
else, x = cheby(n); y = cheby(m); end 
V = ones(n); for j=2:n, V(:,j) = V(:,j-1).*x; end % Vandermonde, original nodes
R = ones(m,n); for j=2:n, R(:,j) = R(:,j-1).*y; end % monomial eval matrix @ y
P = (V'\R')';                                       % backwards-stable solve
end


