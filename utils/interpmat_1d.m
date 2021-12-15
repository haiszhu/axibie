function L = interpmat_1d(t,s)
% INTERPMAT_1D   interpolation matrix from nodes in 1D to target nodes
%
% L = interpmat_1d(t,s) returns interpolation matrix taking values on nodes s
%  to target nodes t. Computed in Helsing style.

% bits taken from qplp/interpmatrix.m from 2013
% Barnett 7/17/16

if nargin==0, test_interpmat_1d; return; end

p = numel(s); q = numel(t); s = s(:); t = t(:);       % all col vecs
n = p; % set the polynomial order we go up to (todo: check why bad if not p)
V = ones(p,n); for j=2:n, V(:,j) = V(:,j-1).*s; end   % polyval matrix on nodes
R = ones(q,n); for j=2:n, R(:,j) = R(:,j-1).*t; end   % polyval matrix on targs
L = (V'\R')'; % backwards-stable way to do it (Helsing) See corners/interpdemo.m
end 