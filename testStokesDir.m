% micro-swimmer test
% 
% arclength reparametrization is not working right now...
% 
% Hai 12/16/20, to be updated with self convergence test 12/17/21

clear all
setup()

flag = 0;

p = 12; np = 32;
if flag
    lam = 0.75; % 0<lam<1
    s.Z = @(t) -(1.5+cos(t)).*(-sin(lam*pi*sin(t)) + 1i*cos(lam*pi*sin(t)));
    s.p = p; s.tpan = linspace(0,pi,np+1)'; 
else
    a = 2; s.p = p; s.np = np; s.tpan = linspace(0,pi,np+1)';
    s.Z = @(t) (a^(-1/3)*sin(pi-t) + 1i*a^(2/3)*cos(pi-t)); % ellipsoid with conserved volume
    s.p = p; s.np = numel(s.tpan)-1; [s,~] = quadr(s, [], 'p', 'G'); 
end
s = quadr(s,[],'p','G');

% target
nx = 100; gx = ((1:nx)/nx*2)*2; ny = 100; gy = ((1:ny)/ny*2-1)*3; % set up plotting grid
[xx, yy] = meshgrid(gx,gy); zz = (xx+1i*yy);
t = [];
[IN, ON] = inpolygon(real(zz),imag(zz),real(s.x),imag(s.x));
ii = (~IN)&(~ON);
t.x = zz(ii(:)); 

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

figure(1),clf,
plot(s.x,'k.-'); axis equal, hold on, quiver(real(zz),imag(zz),real(u),imag(u),'b')


% self test
if 1
    s.tpan = linspace(0,pi,2*np+1)';
    s = quadr(s,[],'p','G');
    f = ([sin(s.t)-1/2*sin(2*s.t);sin(s.t)-1/2*sin(2*s.t)].*[real(s.tang);imag(s.tang)]);
    A = AxiStokesSpecialMat(s,s);
    mu = A\f;
    u2 = nan*(1+1i)*zz;
    At = AxiStokesSpecialMat(t,s);
    temp = At*mu;
    u2(ii) = temp(1:end/2) + 1i*temp(end/2+1:end);
    figure(2),clf,
    imagesc(gx,gy,log10(abs(u-u2))), axis equal
    colorbar, hold on, plot(s.xlo,'or'), hold off,    
    title('log10 err in |u|'),% caxis([-15 -5])
end


keyboard

