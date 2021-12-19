function [A,Az] = S0formChebyquad(t,s,a,b,interpMat,w0f)
% special quadrature based on a differential 0-form approach 
% with cheby basis
%
% Hai 05/11/21

N = numel(t.x);                          % # of targets
p = numel(s.x);                          % assume panel order is # nodes
if nargin<5 
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
idx1 = (abs(xnr+1) + abs(xnr-1)) <= 1.1*2;  % closest
idx2 = ~idx1;    % middle

% closest 
x1 = xnr(idx1); N1 = numel(x1);
T1 = zeros(p+1,N1); 
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
T1(1,:) = p1r'+1i*p1i'; T1(2,:) = x1.'.*T1(1,:) + 2; k=2; T1(k+1,:) = 2*x1.'.*T1(k,:) - T1(k-1,:);
for k=3:p
    T1(k+1,:) = (1/k-1/(k-2))-((-1)^k/k-(-1)^(k-2)/(k-2)) + 2*x1.'.*T1(k,:) - T1(k-1,:);
end

% middle (forward recursion fails, instead a solving process)
x2 = xnr(idx2); N2 = numel(x2);
T2 = zeros(p+1,N2); 
p2i = -atan2(-imag(x2),(real(x2).^2+imag(x2).^2-real(x2)))+atan2(imag(x2),(real(x2).^2+imag(x2).^2+real(x2)));
p2r = 1/2*log(((1-real(x2)).^2+imag(x2).^2)./((1+real(x2)).^2+imag(x2).^2));
T2(1,:) = p2r'+1i*p2i'; T2(2,:) = x2.'.*T2(1,:) + 2; 
yf = interpMat*y; wxp = w0f.*(interpMat*yxp); 
tmp1 = ones(size(yf)); tmp2 = yf; yfp = zeros(size(yf));
for k=3:p, yfp = 2*yf(:).*tmp2-tmp1; tmp1 = tmp2; tmp2 = yfp; end  % Chebyshev 
Tp = (-1./bsxfun(@minus,x2,yf.')*(wxp.*yfp)).'; 
rhs0 = -(1./(3:p-1)-1./(1:p-3)) + ((-1).^(3:p-1)./(3:p-1)-(-1).^(1:p-3)./(1:p-3)); rhs0 = rhs0(:); 
Rhs = rhs0*ones(1,N2); Rhs(1,:) = Rhs(1,:) + T2(2,:); Rhs(end,:) = Rhs(end,:) + Tp; % Rhs = Rhs(:);
h=1/(p-2); E = sin((1:(p-3))'*(1:(p-3))*pi*h); Ei = E\Rhs; Lam = bsxfun(@plus,-2*cos((1:(p-3))'*pi*h),2*x2.'); 
T2(3:p-1,:) = E*(Ei./Lam); T2(p,:) = Tp;
yfp = 2*yf(:).*tmp2-tmp1; Tp = (-1./bsxfun(@minus,x2,yf.')*(wxp.*yfp)).'; T2(p+1,:) = Tp;

% something special about SLP
T = zeros(p+1,N1+N2); T(:,idx1) = T1; T(:,idx2) = T2;
Q = zeros(p,N1+N2); 
Q(1,:) = T(2,:) - log((1-xnr.').*(-1-xnr.'));
Q(2,:) = 1/4*(T(3,:) - T(1,:));
Q(3,:) = 1/2*(1/3*T(4,:)-T(2,:)) - 1/2*(1/3-1)*log((1-xnr.').*(-1-xnr.'));
Q(4,:) = 1/2*(1/4*T(5,:)-1/2*T(3,:)) - 1/2*(1/4-1/2)*T(1,:);
Q(3:2:end,:) = 1/2*(T(4:2:end,:).*repmat(1./(3:2:p)',[1 N1+N2]) - T(2:2:end-2,:).*repmat(1./(1:2:p-2)',[1 N1+N2]))...
                -1/2*(1./(3:2:p)'-1./(1:2:p-2)')*log((1-xnr.').*(-1-xnr.'));
Q(4:2:end,:) = 1/2*(T(5:2:end,:).*repmat(1./(4:2:p)',[1 N1+N2]) - T(3:2:end-2,:).*repmat(1./(2:2:p-2)',[1 N1+N2]))...
                -1/2*(1./(4:2:p)'-1./(2:2:p-2)')*T(1,:);

% matrix
Anr = real((VT.'\Q).'.*repmat((1i*s.nx)',[N1+N2 1])*zsc)/(2*pi*abs(zsc)); 
Anr = Anr*abs(zsc) - log(abs(zsc))/(2*pi)*repmat(abs(s.wxp)',[N1+N2 1]); % unscale, yuk
rfr = bsxfun(@minus,xfr,yf.'); wspfr = bsxfun(@times,abs(wxp),interpMat);
Afr = -(1/2/pi)*log(abs(rfr))*wspfr; % upsampling for far
A = zeros(N,p); A(idxnr,:) = Anr; A(idxfr,:) = Afr;

% derivative matrix
Aznr = (VT.'\T(1:end-1,:)).'*(1/(2*pi)).*repmat((1i*s.nx)',[N1+N2 1]);  % solve for targ complex-deriv mat & rescale
wxpfr = bsxfun(@times,wxp,interpMat);
Azfr = -(1/(2*pi)./rfr*wxpfr).*repmat((1i*s.nx)',[N-N1-N2 1]); % upsampling for far
Az = zeros(N,p); Az(idxnr,:) = Aznr; Az(idxfr,:) = Azfr;

% keyboard

end