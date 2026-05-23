function A = LapDLPAxiMat(t,s)
% Axisymmetric Laplace double-layer (DLP) matrix, naive (smooth) quadrature.
%
% Modal n=0 double layer  D[mu](r) = int_Gamma dG/dn(y) mu(y) dS_y,
%   dG/dn(y) = (1/4pi) n_y.(x-y)/|x-y|^3      (normal at the SOURCE y).
%
% Far/smooth operator only: all kernels treated as smooth (weights s.ws);
% close/self interactions need a special quadrature.  This mirrors
% LapGradAxiMat.m (the S'/adjoint DLP) with the target normal replaced by the
% source normal:  n_rho/rho -> n_rho'/rho',  m -> m', and the Cauchy (E) term
% sign flipped (- -> +).  See latex/laplace_dlp_axi.tex (blocks B04, B09).
%
% s.x source, s.nx source normal, s.ws source weights, t.x targets.
%
% Hai (DLP analog of LapGradAxiMat)

x1 = real(t.x(:)); x2 = imag(t.x(:));        % target (column)
y1 = real(s.x(:))'; y2 = imag(s.x(:))';      % source (row)
ny1 = real(s.nx(:))'; ny2 = imag(s.nx(:))';  % source normal (row)
r2 = bsxfun(@minus,x1,y1).^2+bsxfun(@minus,x2,y2).^2;
chi = 1+r2./(2*bsxfun(@times,x1,y1)); tmp = 2./(chi+1);
[SmatK,SmatE] = ellipke(tmp);
SnaiveK = sqrt(tmp).*SmatK; SnaiveE = sqrt(tmp).*SmatE;
sqrty1x1 = sqrt(bsxfun(@rdivide,y1,x1)); x1y1 = bsxfun(@times,x1,y1);
% s_K = -sqrt(rho'/rho) * n_rho'/rho'       (note rho' in denominator, source normal)
SK = -sqrty1x1.*bsxfun(@times,ny1,x1)./x1y1;
% m' = n_rho'(rho-rho') + n_zeta'(zeta-zeta')   (source-normal version of m)
mp = bsxfun(@minus,bsxfun(@times,ny1,x1),bsxfun(@times,ny1,y1)) ...
   + bsxfun(@minus,bsxfun(@times,ny2,x2),bsxfun(@times,ny2,y2));
% s_E = sqrt(rho'/rho) ( n_rho'/rho' + 2 m'/|r-r'|^2 )
%     = 2 sqrt(rho'/rho) ( n_rho' rho (chi-1) + m' ) / |r-r'|^2     (+ sign vs S')
SE = 2*sqrty1x1.*(bsxfun(@times,ny1,x1).*(chi-1)+mp)./r2;
A = 1/(4*pi)*(SnaiveK.*SK + SnaiveE.*SE)*diag(s.ws);

end
