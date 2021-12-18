function PMat = SphereSLPMatP(s,epsilon,tol)
%
%
% Hai 06/15/20
if nargin < 2, epsilon = 1e-05; end
if nargin < 3, tol = 1e-05; end

PMat = NaN(numel(s.x),2*numel(s.x));
t_in0 = [s.tlo;s.thi(end)]; 


%% loop over targets, and for each target setup adaptive grid...(slow)

t_in = t_in0; 
% refine 1st panel
lam = 1/2; n_split = 1;
while t_in(2) > s.t(1)
    split = (1-lam)*t_in(1) + lam*t_in(2);
    t_in = [t_in(1); split; t_in(2:end)];
    n_split = n_split+1;
end
p = s.p; s_aux1.Z = s.Z; s_aux1.p = p; s_aux1.tpan = t_in; qtype = 'p'; qntype = 'G';
s_aux1 = quadr(s_aux1, [], qtype,qntype);  % setup a shado grid, so that all targets would have nearby panels
for k=1:p    
    t.x = s.x(k)+epsilon*s.nx(k); t.t = s.t(k);
%     t.x = s.x(k); t.t = s.t(k);

    i_aux = 1; while t.t > s_aux1.thi(i_aux), i_aux=i_aux+1; end     % which panel is closest to t.x
    t_in_i = [t_in(1:i_aux-1);t_in(i_aux+2:end)];   % before adaptivity, combine near neighbors
    idx = (i_aux-2)*p+(1:3*p);

    t_close = [t_in(i_aux-1);t.t;t_in(i_aux+2)];
    % refine left panel
    lam = 1/2; n_split = 1;
    while (t.t-t_close(n_split)) > tol
        split = (1-lam)*t_close(n_split) + lam*t.t;
        t_close = [t_close(1:n_split); split; t_close(n_split+1:end)];
        n_split = n_split+1;
    end
    % refine right panel
    n_split = n_split+2;
    while tol < t_close(n_split)-t.t
        split = (1-lam)*t_close(n_split) + lam*t.t;
        t_close = [t_close(1:n_split-1); split; t_close(n_split:end)];
    end

    % close panel quadr
    p = s.p; s_close.Z = s.Z; s_close.p = p; s_close.tpan = unique([t_in_i;t_close]); qtype = 'p'; qntype = 'G';
    s_close = quadr(s_close, [], qtype,qntype);

    G = SphereQuadtpP(s_close,t);

    % interpolation
    tempt = s_close.t(1:end-p*(s.np-2));
    [stdt, ~] = gauss(p); stdbw = baryweights(stdt);
    idx1 = ( tempt <= s.thi(1)); 
    Lp_up1 = baryprojs(s.t(1:p)-1/2*(s.tlo(1)+s.thi(1)), stdbw,tempt(idx1)-1/2*(s.tlo(1)+s.thi(1)));
    idx2 = ( tempt > s.thi(1))&( tempt <= s.thi(2)); 
    Lp_up2 = baryprojs(s.t(p+(1:p))-1/2*(s.tlo(2)+s.thi(2)), stdbw,tempt(idx2)-1/2*(s.tlo(2)+s.thi(2)));
    PMat(k,:) = G*blkdiag(Lp_up1,Lp_up2,eye(numel(s.x)-2*p),Lp_up1,Lp_up2,eye(numel(s.x)-2*p));    

end

for k=(p+1):(numel(s.x)-p)
    t.x = s.x(k)+epsilon*s.nx(k); t.t = s.t(k);
    
    i = 1; while t.t > s.thi(i), i=i+1; end 
    t_in_i = [t_in0(1:i-1);t_in0(i+2:end)];   % before adaptivity, combine near neighbors
    idx = (i-2)*p+(1:3*p);
    
    t_close = [t_in0(i-1);t.t;t_in0(i+2)];
    % refine left panel
    lam = 1/2; n_split = 1;
    while (t.t-t_close(n_split)) > tol
        split = (1-lam)*t_close(n_split) + lam*t.t;
        t_close = [t_close(1:n_split); split; t_close(n_split+1:end)];
        n_split = n_split+1;
    end
    % refine right panel
    n_split = n_split+2;
    while tol < t_close(n_split)-t.t
        split = (1-lam)*t_close(n_split) + lam*t.t;
        t_close = [t_close(1:n_split-1); split; t_close(n_split:end)];
    end
    
    % close panel quadr
    p = s.p; s_close.Z = s.Z; s_close.p = p; s_close.tpan = unique([t_in_i;t_close]); qtype = 'p'; qntype = 'G';
    s_close = quadr(s_close, [], qtype,qntype);
    
    G = SphereQuadtpP(s_close,t);
    
    % interpolation
    tempt = s_close.t((i-2)*p+(1:p*(s_close.np-s.np+3)));
    [stdt, ~] = gauss(p); stdbw = baryweights(stdt);
    idx1 = ( tempt <= s.thi(i-1)); 
    Lp_up1 = baryprojs(s.t((i-2)*p+(1:p))-1/2*(s.tlo(i-1)+s.thi(i-1)), stdbw,tempt(idx1)-1/2*(s.tlo(i-1)+s.thi(i-1)));
    idx2 = ( tempt <= s.thi(i))&( tempt > s.thi(i-1));
    Lp_up2 = baryprojs(s.t((i-1)*p+(1:p))-1/2*(s.tlo(i)+s.thi(i)), stdbw,tempt(idx2)-1/2*(s.tlo(i)+s.thi(i)));
    idx3 = ( tempt <= s.thi(i+1))&( tempt > s.thi(i));
    Lp_up3 = baryprojs(s.t((i)*p+(1:p))-1/2*(s.tlo(i+1)+s.thi(i+1)), stdbw,tempt(idx3)-1/2*(s.tlo(i+1)+s.thi(i+1)));
    PMat(k,:) = G*blkdiag(eye((i-2)*p),Lp_up1,Lp_up2,Lp_up3,eye(numel(s.x)-(i+1)*p),eye((i-2)*p),Lp_up1,Lp_up2,Lp_up3,eye(numel(s.x)-(i+1)*p));
end

% refine last panel
t_in = t_in0; 
while t_in(end-1) < s.t(end)
    split = (1-lam)*t_in(end) + lam*t_in(end-1);
    t_in = [t_in(1:end-1); split; t_in(end)];
end
p = s.p; s_auxl.Z = s.Z; s_auxl.p = p; s_auxl.tpan = t_in; qtype = 'p'; qntype = 'G';
s_auxl = quadr(s_auxl, [], qtype,qntype);  % setup a shado grid, so that all targets would have nearby panels

for k=(numel(s.x)-p+1):numel(s.x)    
    t.x = s.x(k)+epsilon*s.nx(k); t.t = s.t(k);
    
    i_aux = 1; while t.t > s_auxl.thi(i_aux), i_aux=i_aux+1; end     % which panel is closest to t.x
    t_in_i = [t_in(1:i_aux-1);t_in(i_aux+2:end)];   % before adaptivity, combine near neighbors
    idx = (i_aux-2)*p+(1:3*p);
    
    t_close = [t_in(i_aux-1);t.t;t_in(i_aux+2)];
    % refine left panel
    lam = 1/2; n_split = 1;
    while (t.t-t_close(n_split)) > tol
        split = (1-lam)*t_close(n_split) + lam*t.t;
        t_close = [t_close(1:n_split); split; t_close(n_split+1:end)];
        n_split = n_split+1;
    end
    % refine right panel
    n_split = n_split+2;
    while tol < t_close(n_split)-t.t
        split = (1-lam)*t_close(n_split) + lam*t.t;
        t_close = [t_close(1:n_split-1); split; t_close(n_split:end)];
    end

    % close panel quadr
    p = s.p; s_close.Z = s.Z; s_close.p = p; s_close.tpan = unique([t_in_i;t_close]); qtype = 'p'; qntype = 'G';
    s_close = quadr(s_close, [], qtype,qntype);

    G = SphereQuadtpP(s_close,t);
    
    % interpolation
    tempt = s_close.t(numel(s.x)-2*p+1:end);
    [stdt, ~] = gauss(p); stdbw = baryweights(stdt);
    idx1 = ( tempt <= s.thi(end-1)); 
    Lp_up1 = baryprojs(s.t(numel(s.x)-2*p+(1:p))-1/2*(s.tlo(end-1)+s.thi(end-1)), stdbw,tempt(idx1)-1/2*(s.tlo(end-1)+s.thi(end-1)));
    idx2 = ( tempt > s.thi(end-1)); 
    Lp_up2 = baryprojs(s.t(numel(s.x)-p+(1:p))-1/2*(s.tlo(end)+s.thi(end)), stdbw,tempt(idx2)-1/2*(s.tlo(end)+s.thi(end)));
    PMat(k,:) = G*blkdiag(eye(numel(s.x)-2*p),Lp_up1,Lp_up2,eye(numel(s.x)-2*p),Lp_up1,Lp_up2);
    
end

end