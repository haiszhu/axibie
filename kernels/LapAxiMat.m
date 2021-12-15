function G = LapAxiMat(t,s)
% s.x source, t.x targets
%
% 
% Hai 05/18/21

x1 = real(t.x(:)); x2 = imag(t.x(:));    % target (colume)
y1 = real(s.x(:))'; y2 = imag(s.x(:))';   % source (row)
r2 = bsxfun(@minus,x1,y1).^2+bsxfun(@minus,x2,y2).^2;
chi = 1+r2./(2*bsxfun(@times,x1,y1)); tmp = 2./(chi+1);
[SmatK,~] = ellipke(tmp);
SnaiveK = sqrt(tmp).*SmatK; 
sqrty1x1 = sqrt(bsxfun(@rdivide,y1,x1));
S = sqrty1x1*2;
G = 1/(4*pi)*(SnaiveK.*S)*diag(s.ws);

end