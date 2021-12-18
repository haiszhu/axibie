function chunkSLP = ChunkSphereSLPMatP(s)
% with a different singular quadrature...
%
% Hai 06/22/20

t_in = [s.tlo;s.thi(end)]; 

% refine 1st panel
lam = 1/2; n_split = 1;
while t_in(2) > s.t(1)
    split = (1-lam)*t_in(1) + lam*t_in(2);
    t_in = [t_in(1); split; t_in(2:end)];
    n_split = n_split+1;
end
% refine last panel
while t_in(end-1) < s.t(end)
    split = (1-lam)*t_in(end) + lam*t_in(end-1);
    t_in = [t_in(1:end-1); split; t_in(end)];
end
% aux panel quadr
p = s.p; s_aux.Z = s.Z; s_aux.p = p; s_aux.tpan = t_in; qtype = 'p'; qntype = 'G';
s_aux = quadr(s_aux, [], qtype); Np_aux = numel(s_aux.xlo);

% alpert
[stdt, ~] = gauss(p); [xt0, wt0, ~] = gauss(p);
stdbw = baryweights(stdt); % universal Lagrange weights for panel on [-1,1]
tpanlen = s.thi'-s.tlo';
t = zeros(p,s.np); for i=1:s.np, t(:,i) = s.tlo(i) + (stdt+1)/2*tpanlen(i); end
tper = s.thi(end) - s.tlo(1); % total length of the periodic param domain
tInfl = [t(:,end)-tper, t, t(:,1)+tper]; % "inflate" one panel at each end (for bary interp)
tpanlenInfl = [tpanlen(end), tpanlen, tpanlen(1)]; % same for panel lengths

%% middle
A.s = cell(1,s.np); %A.i = cell(1,s.np); A.j = cell(1,s.np);
for i = 2:s.np-1
    ip = i+1; im = i-1; % +,- neighbor indices
    te(1) = s.tlo(im); te(2) = s.tlo(i); te(3) = s.tlo(ip); te(4) = s.thi(ip); % t-endpoints of 3 panels (no wrapping to periodicity!)
    indP = (i-2)*p+1:(i+1)*p; % source indices of 3p-long neighbor
    xP = s.x(indP); % 3p-long quadr nodes, for interp to Alpert nodes
    spP = s.sp(indP); % 3p-long speeds, for interp to Alpert nodes
    xce = 1/2*(te(1:3)+te(2:4));
    
    A.s{i} = zeros(s.p,6*s.p); %A.i{i} = A.s{i}; A.j{i} = A.s{i};
    for j = 1:p
        indj = (i-1)*p+j; % global index of the j-th local node
        if j <= 5
            [xs0,ws0]=hqsuppquad10(5);
            xs0 = (te(3)-te(2))*xs0/2;
            xend = [te(1), te(2)+(s.t(indj)-xce(2)-(te(3)-te(2))*xt0(5)/2), te(3)+(s.t(indj)-xce(2)-(te(3)-te(2))*xt0(5)/2), te(4)];
            xs0 = xs0+(s.t(indj)-(te(3)-te(2))*xt0(5)/2); % new center
        else
            [xs0,ws0]=hqsuppquad10(6);
            xs0 = (te(3)-te(2))*xs0/2; 
            xend = [te(1), te(2)+(s.t(indj)-xce(2)-(te(3)-te(2))*xt0(6)/2), te(3)+(s.t(indj)-xce(2)-(te(3)-te(2))*xt0(6)/2), te(4)];
            xs0 = xs0+(s.t(indj)-(te(3)-te(2))*xt0(6)/2); % new center
        end
        
        ws0 = (te(3)-te(2))*ws0/2;
    
        % loop over chunkmatc panels... and combine 'at' and 'aw'
        at=[]; aw=[];
        for k=1:3 
            xlo = xend(k); xhi = xend(k+1); % this is usually different from s_aux tlo and thi...
            if k==1 || k==3
                xs = (xhi-xlo)*xt0/2+1/2*(xlo+xhi); ws = (xhi-xlo)*wt0/2;
            else
                xs = xs0; ws = ws0;
            end
            at = [at;xs]; aw = [aw;ws];
        end
        
        tA = struct('x',s.x(indj)); % target: j-th local node
        potAx = zeros(1,3*p); potAy = potAx;
        
        for c=1:3 % loop over 3 local panels of sources (interp from same panel)
            ic = i+c-1; % index of this panel in the inflated interp domain tInfl
            iA = at>te(c) & at<te(c+1); % Alpert indices in this src panel
            % Lagrange basis for panel nodes eval at A nodes in this src panel...
            indc = (c-1)*p+1:c*p;
            L = baryprojs(tInfl(:,ic)*2/tpanlenInfl(ic), stdbw, at(iA)*2/tpanlenInfl(ic));
            % geom setup at Alpert nodes
            xc = xP(indc); spc = spP(indc);
            ax = baryeval(tInfl(:,ic)*2/tpanlenInfl(ic), stdbw, xc, at(iA)*2/tpanlenInfl(ic));
            asp = baryeval(tInfl(:,ic)*2/tpanlenInfl(ic), stdbw, spc, at(iA)*2/tpanlenInfl(ic));
            aws = aw(iA).*asp; 
            sA = struct('x',ax,'ws',aws,'w',aw(iA),'sp',asp); % Alpert sourse;
%             matA = SLPmatrix(tA,sA);
            matA = SphereQuadtpP(sA,tA);
            potAx(:,(c-1)*p+1:c*p) = matA(:,1:end/2)*L; 
            potAy(:,(c-1)*p+1:c*p) = matA(:,end/2+1:end)*L;
            
%             keyboard
        end
%         rows=[indj,N+indj]'*ones(1,6*s.p);columns=ones(2,1)*[indP,N+indP];
        A.s{i}(j,:) = [potAx,potAy]; % A.i{i}([j;j+p],:) = rows; A.j{i}([j;j+p],:) = columns;
        
%         keyboard
%         figure(),plot(s.t(indP),0*xP,'.'); hold on, plot(s.t(j+(i-1)*p),0,'*'),plot(xend,0*xend,'o'),plot(xs0,0*xs0,'.')
    end
    
    
end

%% 1st panel
i = 1; 
t1 = struct('x',s.x(1:p)); s_aux1 = struct('x',s_aux.x(1:(n_split+1)*p),'ws',s_aux.ws(1:(n_split+1)*p)); % Alpert source;
% A_aux1 = SLPmatrix(t1,s_aux1);  A_aux1_0 = A_aux1;  
A_aux1 = SphereQuadtpP(s_aux1,t1); A_aux1_0 = A_aux1; 
for j=1:p
    indj = (i-1)*p+j; % global index of the j-th local node
    % figure out which panel it is located in s_aux
    i_aux = 1; while s.t(indj) > s_aux.thi(i_aux), i_aux=i_aux+1; end
    ip_aux = i_aux+1; im_aux = i_aux-1; % +,- neighbor indices
    te(1) = s_aux.tlo(im_aux); te(2) = s_aux.tlo(i_aux); te(3) = s_aux.tlo(ip_aux); te(4) = s_aux.thi(ip_aux); % t-endpoints of 3 panels (no wrapping to periodicity!)
    indP = (i_aux-2)*p+1:(i_aux+1)*p; % source indices of 3p-long neighbor
    xP = s_aux.x(indP); % 3p-long quadr nodes, for interp to aux nodes
    spP = s_aux.sp(indP); % 3p-long speeds, for interp to aux nodes
    xce = 1/2*(te(1:3)+te(2:4));
    
    idx0 = (i_aux-1)*p; Idx = 1:Np_aux*p; idxc = Idx([1:idx0-p,idx0+2*p+1:Np_aux*p]); 
    if s.t(indj) <= 1/2*(s_aux.t((i_aux-1)*p+5)+s_aux.t((i_aux-1)*p+6))
        [xs0,ws0]=hqsuppquad10(5);
        xs0 = (te(3)-te(2))*xs0/2; % scale 1st, shift later(need to match singularity instead of end points or center), won't break if lam=1/2
        xend = [te(1), te(2)+(s.t(indj)-xce(2)-(te(3)-te(2))*xt0(5)/2), te(3)+(s.t(indj)-xce(2)-(te(3)-te(2))*xt0(5)/2), te(4)];
        xs0 = xs0+(s.t(indj)-(te(3)-te(2))*xt0(5)/2); % new center
        
%         keyboard
    else
        [xs0,ws0]=hqsuppquad10(6);
        xs0 = (te(3)-te(2))*xs0/2; 
        xend = [te(1), te(2)+(s.t(indj)-xce(2)-(te(3)-te(2))*xt0(6)/2), te(3)+(s.t(indj)-xce(2)-(te(3)-te(2))*xt0(6)/2), te(4)];
        xs0 = xs0+(s.t(indj)-(te(3)-te(2))*xt0(6)/2); % new center
        
%         keyboard
    end
    ws0 = (te(3)-te(2))*ws0/2;
    
    % loop over chunkmatc panels... and combine 'at' and 'aw'
    at=[]; aw=[];
    for k=1:3 
        xlo = xend(k); xhi = xend(k+1); % this is usually different from s_aux tlo and thi...
        if k==1 || k==3
            xs = (xhi-xlo)*xt0/2+1/2*(xlo+xhi); ws = (xhi-xlo)*wt0/2;
        else
            xs = xs0; ws = ws0;
        end
        at = [at;xs]; aw = [aw;ws];
    end
    
    tA = struct('x',s.x(indj)); % target: j-th local node
    potAx = zeros(1,3*p); potAy = potAx;
    
    % loop over aux panels...
    for c=1:3 
        ic = i_aux+c-2; % index of this panel
        iA = at>te(c) & at<te(c+1); % Alpert indices in this src panel
        % Lagrange basis for panel nodes eval at A nodes in this src panel...
        indc = (c-1)*p+1:c*p;
        L = baryprojs(s_aux.t((ic-1)*p+(1:p))*2/(s_aux.thi(ic)-s_aux.tlo(ic)), stdbw, at(iA)*2/(s_aux.thi(ic)-s_aux.tlo(ic)));
        % geom setup at extra nodes
        xc = xP(indc); spc = spP(indc);
        ax = baryeval(s_aux.t((ic-1)*p+(1:p))*2/(s_aux.thi(ic)-s_aux.tlo(ic)), stdbw, xc, at(iA)*2/(s_aux.thi(ic)-s_aux.tlo(ic)));
        asp = baryeval(s_aux.t((ic-1)*p+(1:p))*2/(s_aux.thi(ic)-s_aux.tlo(ic)), stdbw, spc, at(iA)*2/(s_aux.thi(ic)-s_aux.tlo(ic)));
        aws = aw(iA).*asp; 
        sA = struct('x',ax,'ws',aws,'w',aw(iA),'sp',asp);  %  sourse;
        matA = SphereQuadtpP(sA,tA);
        potAx(:,(c-1)*p+1:c*p) = matA(:,1:end/2)*L; 
        potAy(:,(c-1)*p+1:c*p) = matA(:,end/2+1:end)*L;
    end
    A_aux1(j,[indP,numel(s_aux1.x)+indP]) = [potAx,potAy];
end
Stdbw = repmat(stdbw/2,1,2);  Stdbw = Stdbw(:);   % 1st and 2nd panel from s
Lp_up1 = baryprojs(s.t(1:2*p)*2/(s.thi(2)-s.tlo(1)), Stdbw,s_aux.t(1:(n_split+1)*p)*2/(s.thi(2)-s.tlo(1)));
A_aux1 = A_aux1*blkdiag(Lp_up1,Lp_up1);

%% last panel
i = s.np; % last panel
t2 = struct('x',s.x(end-p+1:end)); s_aux2 = struct('x',s_aux.x(end-(n_split+1)*p+1:end),'ws',s_aux.ws(end-(n_split+1)*p+1:end)); % Alpert source;
% A_aux2 = SLPmatrix(t2,s_aux2);  % A_aux2_0 = A_aux2; 
A_aux2 = SphereQuadtpP(s_aux2,t2); A_aux2_0 = A_aux2; 
for j=1:p
    indj = (i-1)*p+j; % global index of the j-th local node
    % figure out which panel it is located in s_aux
    i_aux = 1; while s.t(indj) > s_aux.thi(i_aux), i_aux=i_aux+1; end
    ip_aux = i_aux+1; im_aux = i_aux-1; % +,- neighbor indices
    te(1) = s_aux.tlo(im_aux); te(2) = s_aux.tlo(i_aux); te(3) = s_aux.tlo(ip_aux); te(4) = s_aux.thi(ip_aux); % t-endpoints of 3 panels (no wrapping to periodicity!)
    indP = (i_aux-2)*p+1:(i_aux+1)*p; % source indices of 3p-long neighbor
    xP = s_aux.x(indP); % 3p-long quadr nodes, for interp to Alpert nodes
    spP = s_aux.sp(indP); % 3p-long speeds, for interp to Alpert nodes
    xce = 1/2*(te(1:3)+te(2:4));
    
    idx0 = (i_aux-1)*p; Idx = 1:Np_aux*p; idxc = Idx([1:idx0-p,idx0+2*p+1:Np_aux*p]); 
    if s.t(indj) <= 1/2*(s_aux.t((i_aux-1)*p+5)+s_aux.t((i_aux-1)*p+6))
        [xs0,ws0]=hqsuppquad10(5);
        xs0 = (te(3)-te(2))*xs0/2; % scale 1st, shift later(need to match singularity instead of end points or center), won't break if lam=1/2
        xend = [te(1), te(2)+(s.t(indj)-xce(2)-(te(3)-te(2))*xt0(5)/2), te(3)+(s.t(indj)-xce(2)-(te(3)-te(2))*xt0(5)/2), te(4)];
        xs0 = xs0+(s.t(indj)-(te(3)-te(2))*xt0(5)/2); % new center
        
%         keyboard
    else
        [xs0,ws0]=hqsuppquad10(6);
        xs0 = (te(3)-te(2))*xs0/2; 
        xend = [te(1), te(2)+(s.t(indj)-xce(2)-(te(3)-te(2))*xt0(6)/2), te(3)+(s.t(indj)-xce(2)-(te(3)-te(2))*xt0(6)/2), te(4)];
        xs0 = xs0+(s.t(indj)-(te(3)-te(2))*xt0(6)/2); % new center
        
%         keyboard
    end
    ws0 = (te(3)-te(2))*ws0/2;
    
    % loop over chunkmatc panels... and combine 'at' and 'aw'
    at=[]; aw=[];
    for k=1:3 
        xlo = xend(k); xhi = xend(k+1); % this is usually different from s_aux tlo and thi...
        if k==1 || k==3
            xs = (xhi-xlo)*xt0/2+1/2*(xlo+xhi); ws = (xhi-xlo)*wt0/2;
        else
            xs = xs0; ws = ws0;
        end
        at = [at;xs]; aw = [aw;ws];
    end
    
    tA = struct('x',s.x(indj)); % target: j-th local node
    potAx = zeros(1,3*p); potAy = potAx;
    
    % loop over aux panels...
    for c=1:3 
        ic = i_aux+c-2; % index of this panel
        iA = at>te(c) & at<te(c+1); % Alpert indices in this src panel
        % Lagrange basis for panel nodes eval at A nodes in this src panel...
        indc = (c-1)*p+1:c*p;
        L = baryprojs(s_aux.t((ic-1)*p+(1:p))*2/(s_aux.thi(ic)-s_aux.tlo(ic)), stdbw, at(iA)*2/(s_aux.thi(ic)-s_aux.tlo(ic)));
        % geom setup at extra nodes
        xc = xP(indc); spc = spP(indc);
        ax = baryeval(s_aux.t((ic-1)*p+(1:p))*2/(s_aux.thi(ic)-s_aux.tlo(ic)), stdbw, xc, at(iA)*2/(s_aux.thi(ic)-s_aux.tlo(ic)));
        asp = baryeval(s_aux.t((ic-1)*p+(1:p))*2/(s_aux.thi(ic)-s_aux.tlo(ic)), stdbw, spc, at(iA)*2/(s_aux.thi(ic)-s_aux.tlo(ic)));
        aws = aw(iA).*asp; 
        sA = struct('x',ax,'ws',aws,'w',aw(iA),'sp',asp);  %  sourse;
        matA = SphereQuadtpP(sA,tA);
        potAx(:,(c-1)*p+1:c*p) = matA(:,1:end/2)*L; 
        potAy(:,(c-1)*p+1:c*p) = matA(:,end/2+1:end)*L;
    end
    A_aux2(j,[indP,numel(s_aux2.x)+indP]-(s.np+n_split-3)*p) = [potAx,potAy];
    
end
Lp_up2 = baryprojs(s.t(end-2*p+1:end)*2/(s.thi(end)-s.tlo(end-1)), Stdbw,s_aux.t(end-(n_split+1)*p+1:end)*2/(s.thi(end)-s.tlo(end-1)));
A_aux2 = A_aux2*blkdiag(Lp_up2,Lp_up2);


chunkSLP = SphereQuadtpP(s,s);
chunkSLP(1:p,[1:2*p,numel(s.x)+(1:2*p)]) = A_aux1;
chunkSLP((s.np-1)*p+(1:p),[(s.np-2)*p+(1:2*p),numel(s.x)+(s.np-2)*p+(1:2*p)]) = A_aux2;

for i = 2:s.np-1
%     keyboard
    chunkSLP((i-1)*p+(1:p),[(i-2)*p+(1:3*p),numel(s.x)+(i-2)*p+(1:3*p)])=A.s{i};
end

% keyboard
% figure(),plot(s_aux.t(indP),0*xP,'.'); hold on, plot(s.t(j),0,'*'),plot(xend,0*xend,'o'),plot(xs0,0*xs0,'.')



end
