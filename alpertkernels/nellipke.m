function [val1, val2] = nellipke(m)

val1 = zeros(size(m)); val2 = val1;
for k = 1:length(m)
    if m(k) < -10^15
        K = 0; E = sqrt(abs(m(k)));
    elseif m(k) < 0
        [K, E] = ellipke(-m(k)./(1 - m(k))); K = 1./sqrt(1-m(k)).*K; E = sqrt(1-m(k)).*E;
    elseif m(k)>=0 && m(k) <= 1
        [K, E] = ellipke(m(k));
    else
       [K, E] = ellipke(1);
    end
    val1(k) = K; val2(k) = E;
end

% ind1 = find(m<0);
% ind2 = find(m>1);
% ind3 = intersect(find(m>=0), find(m<=1)); 
% 
% K = zeros(size(m)); E = K;
% [K(ind1), E(ind1)] = ellipke(-m(ind1)./(1-m(ind1))); K(ind1) = 1./sqrt(1-m(ind1)).*K(ind1); E(ind1) = 1./sqrt(1-m(ind1)).*E(ind1); 
% [K(ind2), E(ind2)] = ellipke(1+0*ind2);
% [K(ind3), E(ind3)] = ellipke(m(ind3));
% 
%         
        