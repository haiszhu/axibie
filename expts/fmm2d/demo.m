%
%  Test Laplace particle FMMs in R^2
%

addpath ..

nsource = 2000;
source = zeros(2,nsource);

idist=2;

if( idist == 1 ),
theta=rand(1,nsource)*pi;
phi=rand(1,nsource)*2*pi;
source(1,:)=.5*cos(phi);
source(2,:)=.5*sin(phi);
end

if( idist == 2 ),
source(1,:)=rand(1,nsource);
source(2,:)=rand(1,nsource);
end

target = source(:,1:nsource);
target(1,:) = target(1,:) + 10;
[ndim,ntarget] = size(target);

% check SLP
charge = rand(1,nsource); dipstr = 0; dipvec = 0*source;
ifcharge=1; ifdipole=0; 
ifpot = 1; ifgrad = 0; ifhess = 0;
ifpottarg = 1; ifgradtarg = 0; ifhesstarg = 0;

iprec=4;
[U]=lfmm2dpart(iprec,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,ifgrad,ifhess,ntarget,target,ifpottarg,ifgradtarg,ifhesstarg);
[F]=l2dpartdirect(nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,ifgrad,ifhess,ntarget,target,ifpottarg,ifgradtarg,ifhesstarg);
if( ifpot ), U.pot=-U.pot/(2*pi); end, if( ifpottarg ), U.pottarg=-U.pottarg/(2*pi); end
if( ifpot ), F.pot=-F.pot/(2*pi); end, if( ifpottarg ), F.pottarg=-F.pottarg/(2*pi); end


d1 = bsxfun(@minus,target(1,:)',source(1,:));   % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,target(2,:)',source(2,:));
logr = log(sqrt(d1.^2+d2.^2));   % r matrix, size M*N
u = -(1/2/pi) * logr;
diff = (u*charge(:)-U.pottarg(:))/max(abs(U.pottarg(:)));
max(abs(diff(:)))

% check DLP (more like 1/r Cauchy type)
s.x = source(1,:)'+1i*source(2,:)'; s.nx = ones(size(s.x))+1i*ones(size(s.x)); s.ws = ones(size(s.x));
iprec=4; sig = rand(1,numel(s.x));
dipstr = sig.*s.nx.';  
ifpot = 1; ifgrad = 0; ifhess = 0;
ifpottarg = 1; ifgradtarg = 0; ifhesstarg = 0;
[U]=zfmm2dpart(iprec,nsource,source,dipstr,ifpot,ifgrad,ifhess,ntarget,target,ifpottarg,ifgradtarg,ifhesstarg);
[F]=z2dpartdirect(nsource,source,dipstr,ifpot,ifgrad,ifhess,ntarget,target,ifpottarg,ifgradtarg,ifhesstarg);

A = LapDLPmat(s,s); A(1:numel(s.x)+1:end) = 0;
uvec = A*sig(:);
diffD = max(abs(real(U.pot(:))/(2*pi)-uvec(:)))

t.x = target(1,:)'+1i*target(2,:)';
At = LapDLPmat(t,s);
uvect = At*sig(:);
diffD2 = max(abs(real(U.pottarg(:))/(2*pi)-uvect(:)))


keyboard