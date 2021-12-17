% test a Laplace Neumann BVP
%
% Hai 05/18/21 updated convergence plot 12/17/21, maybe add more BVP types later

setup()

v = 1;
s.p = 16; 
shape = 'star'; % or 'star', or 'star2'. test convergence on different shapes

% family of shapes...
if strcmp(shape,'ellipse')
    lam = 1; ratio = 2; s.Z = @(t) lam*(sin(t) - ratio*1i*cos(t));
    Np = 4:2:24;
elseif strcmp(shape,'star')
    a = .3; w = 5; R = @(t) 2*(1 + a*cos(w*t-pi/2)); s.Z = @(t) R(t).*exp(1i*t); s.Z = @(t) s.Z(t-pi/2);
    Np = 20:2:40;
elseif strcmp(shape,'star2')
    s.Z = @(t) ((cos(t).^2+9*sin(t).^2).^(1/2)+cos(4*t).^2).*sin(t)...
                -1/2*1i*((cos(t).^2+9*sin(t).^2).^(1/2)+cos(4*t).^2).*cos(t);
    np = 48;
    Np = 26:2:50;
end
tpan = linspace(0,pi,24+1)'; s.tpan = tpan;
s = quadr(s,[],'p','G');

% target
nx = 150; gx = ((1:nx)/nx*2)*4; ny = 150; gy = ((1:ny)/ny*2-1)*6; % set up plotting grid
[xx, yy] = meshgrid(gx,gy); zz = (xx+1i*yy);
t = [];
[IN, ON] = inpolygon(real(zz),imag(zz),real(s.x),imag(s.x));
ii = ~IN;
t.x = zz(ii(:));  % eval pts only on one side
t.nx = ones(size(t.x))/2 + 1i*ones(size(t.x));

% exact soln
alpha = linspace(-pi/3,pi/5,5);
y_source = [];  
y_source.x = 3/5*(1*exp(1i*alpha).'+0.1+1i*0.3); y_source.ws = ones(size(y_source.x));
pt_force = [1;1;-1;1;1/2];
fhom = nan*xx; % exact soln
Ahom = LapAxiMat(t,y_source);
fhom(ii(:)) = Ahom*pt_force;

% h refinement, increase num of panels
err = NaN(size(Np));
for k=1:numel(Np)
    
    % boundary discretization
    np = Np(k);
    tpan = linspace(0,pi,np+1)'; s.tpan = tpan;
    s = quadr(s,[],'p','G');

    % use grad SLP to solve a Neumann problem 
    A = LapGradAxiSpecialMat(s,s); % A = (-1/2*I+S')
    rhs = LapGradAxiMat(s,y_source)*pt_force;
    tau = A\rhs;

    % use SLP to evaluate soln
    u = nan*xx; % soln
    u(ii(:)) = LapAxiSpecialMat(t,s)*tau;

    % plot error
    figure(1),clf,imagesc(gx,gy,log10(abs(u-fhom))), 
    colorbar, hold on, fill(real([s.x;s.x(1)]),imag([s.x;s.x(1)]),'w')
    plot(s.x,'-k'), plot([s.xlo;s.xhi(end)],'.r','MarkerSize',12); caxis([-15 0]), colormap('jet')
    axis equal tight

    err(k) = max(abs(u(:)-fhom(:)));
end
figure(1),clf,imagesc(gx,gy,log10(abs(u-fhom))), 
colorbar, hold on, fill(real([s.x;s.x(1)]),imag([s.x;s.x(1)]),'w')
plot(s.x,'-k'), plot([s.xlo;s.xhi(end)],'.r','MarkerSize',12); colormap('jet')
title('log10 err in |u|'), axis equal tight

figure(2),clf, semilogy(s.p*Np,err,'o-k'), title('BVP conv')
xlabel('$n_p$','interpreter','latex')
% fit_order = [log(Np'),ones(numel(Np),1)]\log(err'); 
% errfit = @(x) err(1)/Np(1)^(-2*s.p+1)*x.^(-2*s.p+1);
% p4=loglog(Np(:),err(:),'p-k','MarkerFaceColor','k');hold on, loglog(Np(:),errfit(Np(:)),'--k')

keyboard