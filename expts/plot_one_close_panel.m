% Show error for single Laplace S/D close panel eval, vs vast oversampling.
% Barnett 4/4/19, modified by Hai 05/22/21

addpath(genpath('../'));
addpath('fmm2d')
clear
k = 0.25;   % curvature of panel
s.Z = @(t) t + 1i*k*t.^2 - 1i*k;    % map from [-1,1] to panel, endpts at +-1
s.Zp = @(t) 1 + 2i*k*t;             % deriv, must match s.Z
s.p = 16;           % Helsing panel order
s.tlo = -1; s.thi = 1; s = quadr(s,s.p,'p','G');    % build one panel
sigma = @(t) log(2+t).*sin(3*t + 1);     % lap dens (S or D), scalar, real

dx = 0.01; g = (-1.3:dx:1.3);   % eval and plot grid
[xx yy] = meshgrid(g,g); nt = numel(xx); tt = [xx(:)';yy(:)'];  % 2 by nt

lptype = 'd';   % 's' SLP or 'd' DLP.
be = 2.0;     % upsampling factor for Helsing eval  (1.0 is not good)

% fill Lap Helsing special quad mat then matvec it...
side = 'i';
ttcmplx = tt(1,:)'+1i*tt(2,:)';  % must be col vec for stokespanelcorm
%A = stokespanelcorm(ttcmplx, s, s.Z(-1), s.Z(1), lptype, side, 'G');  % special
% there is no upsample-then-special-eval, so we have to do by hand...
t = []; t.x = ttcmplx;       % target struct
sf = quadr_panf(s,be,'G');   % upsampled panel
if lptype=='s'
  [A,Az] = S0formChebyquad(t,sf,s.Z(-1),s.Z(1));
elseif lptype=='d'
  [A0,A0z,A0zz] = Dspecialquad(t,sf,s.Z(-1),s.Z(1),side);
  [A,Az,Azz] = D0formChebyquad(t,sf,s.Z(-1),s.Z(1));
end
densvals = sigma(s.t);
Imn = interpmat(s.p,sf.p,'G'); densvalsf = Imn*densvals;   % upsamp dens
u = A*densvalsf;                 % col vec of potentials
u = reshape(real(u),size(xx));   % note real part for DLP case
if lptype=='d'
    du = Az*densvalsf; du = reshape(du,size(xx));
    du0 = A0z*densvalsf; du0 = reshape(du0,size(xx));
    ddu = Azz*densvalsf; ddu = reshape(ddu,size(xx));
    ddu0 = A0zz*densvalsf; ddu0 = reshape(ddu0,size(xx));
end
if lptype=='s'
    du = Az*densvalsf; du = reshape(du,size(xx));
    du0 = A0z*densvalsf; du0 = reshape(du0,size(xx));
end

ne = 5e3;         % Compute Lap "exact" ("e") : oversampled source
[te,we] = lgwt(ne,-1,1);
iprec=5;   % FMM highest precision
source0 = s.Z(te); source = [real(source0),imag(source0)]';   % oversam nodes
charge = sigma(te).*abs(s.Zp(te)).*we;     % oversamp dens * speedweight
FMMprefac = -1/(2*pi); charge = charge*FMMprefac;
dipstr = charge; dipvec = -1i*s.Zp(te); dipvec=dipvec./abs(dipvec);
dipvec = [real(dipvec),imag(dipvec)]';                     % unit normals
U = rfmm2dpart(iprec,ne,source,lptype=='s',charge,lptype=='d',dipstr,dipvec,0,0,0,nt,tt,1,0,0);
ue = reshape(U.pottarg, size(xx));
if lptype=='s'
    Azn = -(1/(2*pi)./bsxfun(@minus,t.x,source0(:).')*spdiags(s.Zp(te).*we,0,ne,ne)).*repmat(conj(1i*(dipvec(1,:)+1i*dipvec(2,:))),[numel(t.x) 1]);
    due = reshape(Azn*sigma(te),size(xx));
end
if lptype=='d'
    Azn = 1i/(2*pi)./bsxfun(@minus,t.x,source0(:).').^2; chargex = sigma(te).*s.Zp(te).*we;
    due = reshape(Azn*chargex,size(xx));
    Azzn = -1i/(2*pi)./bsxfun(@minus,t.x,source0(:).').^3;
    ddue = reshape(Azzn*chargex,size(xx));
end
if lptype=='s'
    figure; set(gcf,'position',[200 200 1500 500]);
    subplot(1,4,1); imagesc(g,g,ue); axis equal xy tight;
    sc = max(abs(caxis)); caxis(sc*[-1 1]); colorbar;
    hold on; plot(source(1,:),source(2,:),'k-');
    title(sprintf('ue (nearly exact), lp=%c',lptype));
    subplot(1,4,2); imagesc(g,g,u); axis equal xy tight; caxis(sc*[-1 1]);colorbar;
    hold on; plot(source(1,:),source(2,:),'k-');
    title(sprintf('u Helsing side=%c',side));
    subplot(1,4,3); imagesc(g,g,log10(abs(ue-u)));
    axis equal xy tight; caxis([-15 0]); colorbar;
    hold on; plot(source(1,:),source(2,:),'k-'); title('log_{10} err in Helsing u');
    subplot(1,4,4); imagesc(g,g,log10(abs(du-due)));
    axis equal xy tight; caxis([-15 0]); colorbar;
    hold on; plot(source(1,:),source(2,:),'k-'); title('log_{10} err in Helsing \nabla u');
end
if lptype=='d'
    figure; set(gcf,'position',[200 200 1500 500]);
    subplot(1,5,1); imagesc(g,g,ue); axis equal xy tight;
    sc = max(abs(caxis)); caxis(sc*[-1 1]); colorbar;
    hold on; plot(source(1,:),source(2,:),'k-');
    title(sprintf('ue (nearly exact), lp=%c',lptype));
    subplot(1,5,2); imagesc(g,g,u); axis equal xy tight; caxis(sc*[-1 1]);colorbar;
    hold on; plot(source(1,:),source(2,:),'k-');
    title(sprintf('u Helsing side=%c',side));
    subplot(1,5,3); imagesc(g,g,log10(abs(ue-u)));
    axis equal xy tight; caxis([-15 0]); colorbar;
    hold on; plot(source(1,:),source(2,:),'k-'); title('log_{10} err in Helsing u');
    subplot(1,5,4); imagesc(g,g,log10(abs(du-due)));
    axis equal xy tight; caxis([-15 0]); colorbar;
    hold on; plot(source(1,:),source(2,:),'k-'); title('log_{10} err in Helsing \nabla u');
    subplot(1,5,5); imagesc(g,g,log10(abs(ddu-ddue)));
    axis equal xy tight; caxis([-15 0]); colorbar;
    hold on; plot(source(1,:),source(2,:),'k-'); title('log_{10} err in Helsing \nabla^2 u');
end

keyboard

function sf = quadr_panf(s, be, qntype)  
% set up quadrature on a closed segment
% QUADR_panf - set up quadrature (either coarse or fine nodes) on a segment struct
%
% sf = quadr_panf(s, be) gives quadrature on coarse or fine nodes.
% Inputs: s  = segment struct containing parametrization
%         be = factor by which to increase panel nodes
%  
% Outputs: sf - segment struct on fine node
if be == 1
    sf = s;
    if ~isfield(s,'p')
    s.p=16; 
    end
    p = s.p; % default panel order
    if qntype=='G', [~, w, D] = gauss(p); else, [~, w, D] = cheby(p); end 
    sf.xp = D*sf.x;
    sf.xpp = D*sf.xp;   % acceleration Z''(sf.x)
    sf.w = w;
    sf.sp = abs(sf.xp); sf.tang = sf.xp./sf.sp; sf.nx = -1i*sf.tang;    % outward unit norma
else
if ~isfield(s,'p')
    s.p=16; 
end
p = s.p; % default panel order
sf=[]; sf.p=ceil(be*s.p); pf=sf.p;
Imn = interpmat(p, pf, qntype);
sf.x = Imn*s.x;
if qntype=='G', [xx, w, D] = gauss(pf); else, [xx, w, D] = cheby(pf); end 
if ~isfield(s,'Zp') 
    if ~isfield(s,'xp'), sf.xp = D*sf.x;  else, sf.xp = Imn*s.xp*(s.thi-s.tlo)/2; end  % velocities Z'(sf.x)
else
    sf.xp = 1/2*(s.thi-s.tlo)*s.Zp(s.tlo + (1+xx)/2*(s.thi-s.tlo));
end
if ~isfield(s,'Zpp')
    if ~isfield(s,'xpp'), sf.xpp = D*sf.xp;  else, sf.xpp = Imn*s.xpp*(s.thi-s.tlo)/2; end  % acceleration Z''(sf.x)
else
    sf.xpp = 1/2*(s.thi-s.tlo)*s.Zpp(s.tlo + (1+xx)/2*(s.thi-s.tlo));
end
sf.w = w;
sf.sp = abs(sf.xp); sf.tang = sf.xp./sf.sp; sf.nx = -1i*sf.tang;    % outward unit normals
sf.cur = -real(conj(sf.xpp).*sf.nx)./sf.sp.^2;
sf.ws = sf.w.*sf.sp; % speed weights
sf.wxp = sf.w.*sf.xp; % complex speed weights (Helsing's wzp)
end
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

function [x,w]=lgwt(N,a,b)
% LGWT  Legendre--Gauss weights and nodes.
%
% [x,w]=lgwt(N,a,b)
%
% This script is for computing definite integrals using Legendre-Gauss 
% Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
% [a,b] with truncation order N
%
% Suppose you have a continuous function f(x) which is defined on [a,b]
% which you can evaluate at any x in [a,b]. Simply evaluate it at all of
% the values contained in the x vector to obtain a vector f. Then compute
% the definite integral using sum(f.*w);
%
% Outputs appear to be column vectors.

% Written by Greg von Winckel - 02/25/2004
N=N-1;
N1=N+1; N2=N+2;

xu=linspace(-1,1,N1)';

% Initial guess
y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);

% Legendre-Gauss Vandermonde Matrix
L=zeros(N1,N2);

% Derivative of LGVM
Lp=zeros(N1,N2);

% Compute the zeros of the N+1 Legendre Polynomial
% using the recursion relation and the Newton-Raphson method

y0=2;

% Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0))>eps
    
    
    L(:,1)=1;
    Lp(:,1)=0;
    
    L(:,2)=y;
    Lp(:,2)=1;
    
    for k=2:N1
        L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
 
    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);   
    
    y0=y;
    y=y0-L(:,N2)./Lp;
    
end

% Linear map from[-1,1] to [a,b]
x=(a*(1-y)+b*(1+y))/2;      

% Compute the weights
w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;

end

