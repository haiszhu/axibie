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
