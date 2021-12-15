% test a Laplace Neumann BVP
%
% Hai 05/18/21

setup()

v = 1;
s.p = 16; 

if 1
    lam = 1; ratio = 2; s.Z = @(t) lam*(sin(t) - ratio*1i*cos(t));
    np = 16;
else
    s.Z = @(t) ((cos(t).^2+9*sin(t).^2).^(1/2)+cos(4*t).^2).*sin(t)...
                -1/2*1i*((cos(t).^2+9*sin(t).^2).^(1/2)+cos(4*t).^2).*cos(t);
    np = 48;
end

tpan = linspace(0,pi,np+1)'; s.tpan = tpan;
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


% us grad SLP to solve a Neumann problem 
A = LapGradAxiSpecialMat(s,s); % A = (-1/2*I+S')
rhs = LapGradAxiMat(s,y_source)*pt_force;
tau = A\rhs;

% use SLP to evaluate soln
u = nan*xx; % soln
u(ii(:)) = LapAxiSpecialMat(t,s)*tau;

Err = abs(u-fhom);
figure(1),clf,imagesc(gx,gy,log10(Err)), 
axis equal; hold on; plot(s.xlo,'o'); colorbar


keyboard

% s.nx = ones(size(s.nx))+1i*zeros(size(s.nx));
% Arho = LapGradAxiSpecialMat(s,s); 
% s.nx = zeros(size(s.nx))+1i*ones(size(s.nx));
% Az = LapGradAxiSpecialMat(s,s);
% s = quadrp(s,[],'p','G');
% A = bsxfun(@times,real(s.nx(:)),Arho) + bsxfun(@times,imag(s.nx(:)),Az);
% 
% trho.x = s.x+1e-5;
% u1 = LapAxiSpecialMat(s,s)*tau;
% u2 = LapAxiSpecialMat(trho,s)*tau;
% test1 = (u2-u1)/1e-04;
% test2 = Arho*tau;
% diff = test1-test2;

