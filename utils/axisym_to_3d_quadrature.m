function s = axisym_to_3d_quadrature(varargin)
% AXISYM_TO_3D_QUADRATURE  Rotate an axisymmetric meridian curve around the z-axis to
% produce a full 3D smooth-Nystrom surface quadrature struct (s.x, s.w, optional s.nx)
% ready to feed into Sto3dSLPmat, Lap3dSLPmat, or any brute-force 3D summation routine.
%
% Two calling forms:
%
%   s = axisym_to_3d_quadrature(chnkr, n_angles)
%       chnkr     : chunkie meridian chunker (chnkr.r row 1 = rho, row 2 = z;
%                   chnkr.n is meridian unit normal; chnkr.wts is meridian arclength weights).
%       n_angles  : number of uniform azimuthal samples theta_j = 2*pi*(j-1)/n_angles.
%       Returns s with s.x, s.w, s.nx (axisymmetric meridian normal rotated).
%
%   s = axisym_to_3d_quadrature(rho, zc, wl, n_angles)
%       rho, zc   : Ngen x 1 generating-curve coordinates (cylindrical rho and z).
%       wl        : Ngen x 1 meridian arclength quadrature weights.
%       n_angles  : number of uniform azimuthal samples.
%       Returns s with s.x, s.w (no s.nx; analytic generating curve carries no normals).
%
% In both forms s.x is 3 x (Nmer * n_angles), s.w is 1 x (Nmer * n_angles), ordered
% meridian-node fastest then azimuth slowest.  The surface element is
%   dS = (meridian arclength weight) * rho * (2*pi/n_angles)
% absorbing the cylindrical Jacobian rho and the uniform azimuthal weight.

if nargin == 2
    % ---- chnkr form ----
    chnkr = varargin{1}; n_angles = varargin{2};
    theta = 2*pi*(0:n_angles-1)/n_angles;
    x_3d  = tensorprod(squeeze(chnkr.r(1,:,:)), cos(theta'));
    y_3d  = tensorprod(squeeze(chnkr.r(1,:,:)), sin(theta'));
    z_3d  = tensorprod(squeeze(chnkr.r(2,:,:)), ones(n_angles,1));
    nx_3d = tensorprod(squeeze(chnkr.n(1,:,:)), cos(theta'));
    ny_3d = tensorprod(squeeze(chnkr.n(1,:,:)), sin(theta'));
    nz_3d = tensorprod(squeeze(chnkr.n(2,:,:)), ones(n_angles,1));
    ws_3d = tensorprod(chnkr.wts .* squeeze(chnkr.r(1,:,:)), 2*pi/n_angles*ones(n_angles,1));
    s.x   = [x_3d(:) y_3d(:) z_3d(:)]';
    s.w   = ws_3d(:)';
    s.nx  = [nx_3d(:) ny_3d(:) nz_3d(:)]';
elseif nargin == 4
    % ---- analytic form (rho, zc, wl, n_angles) ----
    rho = varargin{1}(:);  zc = varargin{2}(:);  wl = varargin{3}(:);  n_angles = varargin{4};
    theta = 2*pi*(0:n_angles-1)/n_angles;
    X = rho * cos(theta);  Y = rho * sin(theta);  Z = repmat(zc, 1, n_angles);
    W = (wl .* rho) * (2*pi/n_angles) * ones(1, n_angles);
    s.x = [X(:).'; Y(:).'; Z(:).'];
    s.w = W(:).';
else
    error('axisym_to_3d_quadrature: pass (chnkr, n_angles) or (rho, zc, wl, n_angles).');
end
end
