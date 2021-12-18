function Ip = AxiKernelP(y1, y2, x1, x2)

n = length(y1);
Ip = zeros(n,2);

if (abs(x1) < 10^-14)
    % panel quadrature probably doesn't need this special case...
    keyboard
end


p = 16; N = 32*p;
tol = 1e-07;
[tt,ww] = weight_setup(p,N, tol);


for k = 1:n
    a = y1(k); b = x1; c = y2(k) - x2;
    
%     if sqrt((a-b)^2+c^2) < 1e-02
    if 0
        f1 = @(t) (a-b*cos(t))./(b^2+a^2-2*a*b*cos(t)+c^2).^(3/2);
        f2 = @(t) c./(b^2+a^2-2*a*b*cos(t)+c^2).^(3/2);
        Ip(k,1) = sum(f1(tt).*ww); Ip(k,2) = sum(f2(tt).*ww);
        
    else
        P = sqrt(a^2+b^2+2*a*b+c^2); M = sqrt(a^2+b^2-2*a*b+c^2);
        [K, E] = nellipke(4*a*b/P^2);
        Ip(k,1) = -2/a/P/M^2*(-M^2*K + (b^2 +c ^2 - a^2)*E);
        Ip(k,2) = 4*c*E/P/M^2;
%         Ip(k,1) = -2*c/a/P/M^2*(-M^2*K + (b^2 +c ^2 - a^2)*E);
%         Ip(k,2) = 4*c^2*E/P/M^2;
        
    end

end

Ip = repmat(y1,1,2).*Ip;
