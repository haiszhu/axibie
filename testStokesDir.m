% micro-swimmer test
% 
% arclength reparametrization is not working right now...
% 
% Hai 12/16/20, to be updated with self convergence test 12/17/21

clear all
setup()

shape = 'cshape';

p = 16; s.p = p; 
if strcmp(shape,'ellipse')
    a = 2; 
    s.Z = @(t) (a^(-1/3)*sin(pi-t) + 1i*a^(2/3)*cos(pi-t)); % ellipsoid with conserved volume
    Np = 2:2:24;
elseif strcmp(shape,'cshape')
    lam = 0.75; % 0<lam<1
    s.Z = @(t) -(1.5+cos(t)).*(-sin(lam*pi*sin(t)) + 1i*cos(lam*pi*sin(t)));
    Np = 6:2:24;
elseif strcmp(shape,'star')
    a = .3; w = 5; R = @(t) 2*(1 + a*cos(w*t-pi/2)); s.Z = @(t) R(t).*exp(1i*t); s.Z = @(t) s.Z(t-pi/2);
    Np = 26:2:44;
elseif strcmp(shape,'star2')
    s.Z = @(t) ((cos(t).^2+9*sin(t).^2).^(1/2)+cos(4*t).^2).*sin(t)...
                -1/2*1i*((cos(t).^2+9*sin(t).^2).^(1/2)+cos(4*t).^2).*cos(t);
    np = 48;
    Np = 26:2:52;
end
s.tpan = linspace(0,pi,24+1)'; s = quadr(s,[],'p','G');

% target
nx = 100; gx = ((1:nx)/nx*2)*2; ny = 100; gy = ((1:ny)/ny*2-1)*3; % set up plotting grid
[xx, yy] = meshgrid(gx,gy); zz = (xx+1i*yy);
t = [];
[IN, ON] = inpolygon(real(zz),imag(zz),real(s.x),imag(s.x));
ii = (~IN)&(~ON);
t.x = zz(ii(:)); 

% reference solution on a fine discretization for self test
s.tpan = linspace(0,pi,(Np(end)+10)+1)';
sf = quadr(s,[],'p','G');
Af = AxiStokesSpecialMat(sf,sf);
muf = Af\([sin(sf.t)-1/2*sin(2*sf.t);sin(sf.t)-1/2*sin(2*sf.t)].*[real(sf.tang);imag(sf.tang)]);
uf = nan*(1+1i)*zz;
Atf = AxiStokesSpecialMat(t,sf);
tempf = Atf*muf;
uf(ii) = tempf(1:end/2) + 1i*tempf(end/2+1:end);


% h refinement, increase num of panels
err = NaN(size(Np));
for k=1:numel(Np)

    % boundary discretization
    np = Np(k);
    tpan = linspace(0,pi,np+1)'; s.tpan = tpan;
    s = quadr(s,[],'p','G');

    % slip velocity
    f = ([sin(s.t)-1/2*sin(2*s.t);sin(s.t)-1/2*sin(2*s.t)].*[real(s.tang);imag(s.tang)]);
    
    % forward solve matrix
    A = AxiStokesSpecialMat(s,s);
    
    % solve
    mu = A\f;
    
    % velocity field
    u = nan*(1+1i)*zz;
    At = AxiStokesSpecialMat(t,s);
    temp = At*mu;
    u(ii) = temp(1:end/2) + 1i*temp(end/2+1:end);

    figure(1),clf,imagesc(gx,gy,log10(abs(u-uf))), 
    colorbar, hold on, fill(real([s.x;s.x(1)]),imag([s.x;s.x(1)]),'w')
    plot(s.x,'-k'), plot([s.xlo;s.xhi(end)],'.r','MarkerSize',12); caxis([-15 0]), colormap('jet')
    axis equal tight
    
    err(k) = max(abs(u(:)-uf(:)));
end

figure(1),clf,imagesc(gx,gy,log10(abs(u-uf))), 
colorbar, hold on, fill(real([s.x;s.x(1)]),imag([s.x;s.x(1)]),'w')
plot(s.x,'-k'), plot([s.xlo;s.xhi(end)],'.r','MarkerSize',12); colormap('jet')
title('log10 err in |u|'), axis equal tight

figure(2),clf, semilogy(s.p*Np,err,'o-k'), title('BVP conv')
xlabel('$n_p$','interpreter','latex')



keyboard

