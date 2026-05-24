function A = LapDLPGradAxiMat(t,s)

x1=real(t.x(:)); x2=imag(t.x(:)); nx1=real(t.nx(:)); nx2=imag(t.nx(:));      % target + target normal (col)
y1=real(s.x(:))'; y2=imag(s.x(:))'; ny1=real(s.nx(:))'; ny2=imag(s.nx(:))';  % source + source normal (row)

dxt=bsxfun(@minus,x1,y1);                 % rho - rho'
dzt=bsxfun(@minus,x2,y2);                 % zeta - zeta'
r2 =dxt.^2+dzt.^2;                        % |r-r'|^2
rr =bsxfun(@times,x1,y1);                 % rho rho'
chi=1+r2./(2*rr); tmp=2./(chi+1);
[SmatK,SmatE]=ellipke(tmp);
SnaiveK=sqrt(tmp).*SmatK; SnaiveE=sqrt(tmp).*SmatE;   % sqrt(2/(chi+1)) K, E
sqrty1x1=sqrt(bsxfun(@rdivide,y1,x1));    % sqrt(rho'/rho)

% normal-displacement products  n.d = n_rho(rho-rho')+n_zeta(zeta-zeta'),  nu.d = ...
md =bsxfun(@times,nx1,dxt)+bsxfun(@times,nx2,dzt);    % n . (r-r')   (target normal)
mpd=bsxfun(@times,ny1,dxt)+bsxfun(@times,ny2,dzt);    % nu . (r-r')  (source normal)

% chi-extended product Phi = (n.d - n_rho |r-r'|^2/2rho)(nu.d + n_rho' |r-r'|^2/2rho')
fac1=md -bsxfun(@rdivide,bsxfun(@times,nx1,r2),2*x1);
fac2=mpd+bsxfun(@rdivide,bsxfun(@times,ny1,r2),2*y1);
Phi =fac1.*fac2;

nn=bsxfun(@times,nx1,ny1);                % n_rho n_rho'
mm=bsxfun(@times,nx2,ny2);                % n_zeta n_zeta'
% s_E^(c) numerator factor  g = n_zeta n_rho' rho - n_rho n_zeta' rho' - n_rho n_rho'(zeta-zeta')
g = bsxfun(@times,bsxfun(@times,nx2,ny1),x1) ...
  - bsxfun(@times,bsxfun(@times,nx1,ny2),y1) ...
  - nn.*dzt;

sK = nn./(2*rr) + Phi./(rr.*(chi+1).*r2);
sE = 2*(chi.*nn+mm)./r2 - 8*chi.*Phi./((chi+1).*r2.^2) + 3*dzt.*g./(rr.*r2);

A = 1/(4*pi)*( SnaiveK.*(sqrty1x1.*sK) + SnaiveE.*(sqrty1x1.*sE) )*diag(s.ws);

end
