function [u,ut]=LapSLPFMM2D(s,if_s,t,if_t,sig,iprec)
% Laplace 2d fmm evaluation
%
% Hai 01/20/21

if nargin < 6, iprec=4; end

nsource=length(s.x); ntarget=length(t.x); % s.x, t.x are all complex numbers.
source=[real(s.x(:)),imag(s.x(:))]'; target=[real(t.x(:)),imag(t.x(:))]';
charge=(s.ws(:).*sig(:))'; dipstr = 0; dipvec = 0*source;

ifcharge=1; 
if if_s==1, ifpot=1; else, ifpot=0; end
if if_t==1, ifpottarg=1; else, ifpottarg=0; end
ifdipole=0; ifgrad=0; ifhess=0; ifgradtarg=0; ifhesstarg=0;
if( ntarget == 0 ), ifpottarg = 0; end
[U]=lfmm2dpart(iprec,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,ifgrad,ifhess,ntarget,target,ifpottarg,ifgradtarg,ifhesstarg);
u = []; ut = [];
if( ifpot ), u=-real(U.pot/(2*pi)); end
if( ifpottarg ), ut=-real(U.pottarg/(2*pi)); end


end