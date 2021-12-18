function G = SphereQuadtp(s,t)
X = [real(s.x);imag(s.x)]; m = length(X)/2;
x1 = X(1:m); x2 = X(1+m:2*m);
% w = s.w;

mt = numel(t.x);
G1 = zeros(mt,m); G2 = G1; G3 = G1; G4 = G1;

for ne = 1:mt    % target
    
%     wt = w;
    wt = ones(size(s.ws));
    Ker = AxiKernel(x1, x2, real(t.x(ne)), imag(t.x(ne))); 
    
    G1(ne, :) = (wt.*(Ker(:,2) + Ker(:,3)))'; 
    G2(ne, :) = (wt.*Ker(:,4))';
    
    G3(ne, :) = (wt.*Ker(:,5))';
    G4(ne, :) = (wt.*(Ker(:,1) + Ker(:,6)))';
end
W = 0.5*s.ws/8/pi; 

G = [G1 G2;G3 G4]; 
G = G.*[repmat(W', 2*mt, 1), repmat(W', 2*mt, 1)];
G = 2*G;

end