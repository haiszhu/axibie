function A = LapGradAxiMat(t,s)
% s.x source, t.x targets, t.nx targets normal
%
%
% Hai 05/17/21

x1 = real(t.x(:)); x2 = imag(t.x(:)); nx1 = real(t.nx(:)); nx2 = imag(t.nx(:));    % target (colume)
y1 = real(s.x(:))'; y2 = imag(s.x(:))';   % source (row)
r2 = bsxfun(@minus,x1,y1).^2+bsxfun(@minus,x2,y2).^2;
chi = 1+r2./(2*bsxfun(@times,x1,y1)); tmp = 2./(chi+1);
[SmatK,SmatE] = ellipke(tmp);
SnaiveK = sqrt(tmp).*SmatK; SnaiveE = sqrt(tmp).*SmatE; 
sqrty1x1 = sqrt(bsxfun(@rdivide,y1,x1)); x1y1 = bsxfun(@times,x1,y1);
SK = -sqrty1x1.*bsxfun(@times,nx1,y1)./x1y1; 
SE = 2*sqrty1x1.*(bsxfun(@times,nx1,y1).*(chi-1)-(bsxfun(@minus,nx1.*x1,bsxfun(@times,nx1,y1))+bsxfun(@minus,nx2.*x2,bsxfun(@times,nx2,y2))))./r2;
A = 1/(4*pi)*(SnaiveK.*SK + SnaiveE.*SE)*diag(s.ws);

end