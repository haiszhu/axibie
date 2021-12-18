function I = AxiKernel(y1, y2, x1, x2)

n = length(y1);
I = zeros(n,6);
flag = 1;

if (abs(x1) < 10^-14)
    for k = 1:n
        a = y1(k); b = x1; c = y2(k) - x2;
        P = sqrt(a^2+b^2+2*a*b+c^2); M = sqrt(a^2+b^2-2*a*b+c^2);

        [K, E] = nellipke(0);

        I(k,1) = 4*a*K/P;
        I(k,2) = 0;
        I(k,3) = 0;
        I(k,4) = 0;
        I(k,5) = -2*c/P/M^2*(-M^2*K + (b^2 +c ^2 - a^2)*E);
        I(k,6) = 4*c^2*a*E/P/M^2;
    end
return;
end

if ~flag
    p = 16; N = 32*p;
    tol = 1e-07;
    [tt,ww] = weight_setup(p,N, tol);
end
    
for k = 1:n
    a = y1(k); b = x1; c = y2(k) - x2;
    if flag
        P = sqrt(a^2+b^2+2*a*b+c^2); M = sqrt(a^2+b^2-2*a*b+c^2);

        [K, E] = nellipke(4*a*b/P^2);

        I(k,1) = 4*K/P;
        I(k,2) = 2/a/b/P*( (a^2+b^2+c^2)*K - P^2*E);
        I(k,3) = 2*c^2/a/b/P/M^2*(M^2*K - (a^2+b^2+c^2)*E);
        I(k,4) = 2*c/b/P/M^2*(-M^2*K + (a^2 +c ^2 - b^2)*E);
        I(k,5) = -2*c/a/P/M^2*(-M^2*K + (b^2 +c ^2 - a^2)*E);
        I(k,6) = 4*c^2*E/P/M^2;
    else
        f1 = @(t) 1./(b^2+a^2-2*a*b*cos(t)+c^2).^(1/2);
        f2 = @(t) cos(t)./(b^2+a^2-2*a*b*cos(t)+c^2).^(1/2);
        f3 = @(t) (a*cos(t)-b).*(a-b*cos(t))./(a^2+b^2-2*a*b*cos(t)+c^2).^(3/2);
        f4 = @(t) c*(a*cos(t)-b)./(b^2+a^2-2*a*b*cos(t)+c^2).^(3/2);
        f5 = @(t) c*(a-b*cos(t))./(b^2+a^2-2*a*b*cos(t)+c^2).^(3/2);
        f6 = @(t) c^2./(b^2+a^2-2*a*b*cos(t)+c^2).^(3/2);


        I(k,1) = sum(f1(tt).*ww); I(k,2) = sum(f2(tt).*ww); I(k,3) = sum(f3(tt).*ww);
        I(k,4) = sum(f4(tt).*ww); I(k,5) = sum(f5(tt).*ww); I(k,6) = sum(f6(tt).*ww);
    end
end
I = repmat(y1,1,6).*I;

