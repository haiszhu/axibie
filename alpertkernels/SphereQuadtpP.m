function G = SphereQuadtpP(s,t)
% traction

X = [real(s.x);imag(s.x)]; m = length(X)/2;
x1 = X(1:m); x2 = X(1+m:2*m);
% w = s.w;

mt = numel(t.x);
G1 = zeros(mt,m); G2 = G1; 

for ne = 1:mt    % target
    
%     wt = w;
    wt = ones(size(s.ws));
    Ker = AxiKernelP(x1, x2, real(t.x(ne)), imag(t.x(ne))); 
    
    G1(ne, :) = (wt.*(Ker(:,1)))'; 
    G2(ne, :) = (wt.*Ker(:,2))';
    
end
W = s.ws/4/pi; 

G = [G1 G2]; 
G = G.*[repmat(W', mt, 1), repmat(W', mt, 1)];

end