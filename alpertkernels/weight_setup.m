function [tt,ww] = weight_setup(p,N, tol)

% tol = 1e-09;

np = ceil(N/p);     
tlo = (0:np-1)'/np*2*pi; thi = (1:np)'/np*2*pi; 
t_in = [tlo;thi(end)]; 

% refine 1st panel
lam = 1/2; n_split = 1;
while t_in(2) > tol
    split = (1-lam)*t_in(1) + lam*t_in(2);
    t_in = [t_in(1); split; t_in(2:end)];
    n_split = n_split+1;
end
% refine last panel
while tol < t_in(end)-t_in(end-1)
    split = (1-lam)*t_in(end) + lam*t_in(end-1);
    t_in = [t_in(1:end-1); split; t_in(end)];
end

tlo = t_in(1:end-1); thi = t_in(2:end);
% keyboard

pt = thi - tlo; np = numel(tlo);
tt = zeros(p*np,1); ww = tt; [x, w0] = gauss(p); 
for i=1:np
    ii = (i-1)*p+(1:p); 
    tt(ii) = tlo(i) + (1+x)/2*pt(i); ww(ii) = w0*pt(i)/2; 
end

end