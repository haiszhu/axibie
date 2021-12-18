% test both Dirichlet and Neumann BVP for axi symmetric solver
%
% Hai 04/25/20

setup()

v = 1;
lptype = 's'; % test SLP (Dirichlet) or SLPt (Neumann)

% generating curve 
% (make sure surface normal direction for Neumann or traction computation)
if 0
    lam = 0.95; % 0<lam<1
    s0.Z = @(t) -(1.5+cos(t)).*(-sin(lam*pi*sin(t)) + 1i*cos(lam*pi*sin(t)))+1i*0.85;
elseif 0
    ratio = 1; % if change ratio, need to change stokeslet position accordingly
    s0.Z = @(t) 1.5*(sin(t) - 1i*ratio*cos(t)); 
else 
	s0.Z = @(t) ((cos(t).^2+9*sin(t).^2).^(1/2)+cos(4*t).^2).*sin(t)...
                -1/2*1i*((cos(t).^2+9*sin(t).^2).^(1/2)+cos(4*t).^2).*cos(t);
end

% target
p=10; N = 10*p; qtype = 'p'; qntype = 'G'; s0.p = p;
[s,~] = quadr(s0, N, qtype, qntype); 
nx = 100; gx = ((1:nx)/nx*2)*4; ny = 100; gy = ((1:ny)/ny*2-1)*6; % set up plotting grid
[xx, yy] = meshgrid(gx,gy); zz = (xx+1i*yy);
t = [];
[IN, ON] = inpolygon(real(zz),imag(zz),real(s.x),imag(s.x));
ii = ~IN;
t.x = zz(ii(:));  % eval pts only on one side

% exact soln
alpha = linspace(-pi/3,pi/5,5);
y_force = [];  
if 1
    y_force.x = 2/3*(1*exp(1i*alpha).'+0.1+1i*0.3); 
else
    y_force.x = 0.4*exp(1i*alpha).'+0.1+1i*0.3; 
end
pt_force = [[1;1;0;1;0];[1;0;-1;0;1]];
fhom = nan*(1+1i)*zz; % exact soln
Ahom = AxiStokeslet(y_force,t);
f_temp = Ahom*pt_force;
fhom(ii(:)) = f_temp(1:end/2) + 1i*f_temp(end/2+1:end); 

fhompres = nan*zz; % exact soln
Ahompres = AxiStokesletP(y_force,t);
fpres_temp = Ahompres*pt_force;
fhompres(ii(:)) = fpres_temp; 


% plot the exact soluntion and the point forces
if v==1
    figure(3),
    clf,streamslice(gx,gy,real(fhom),imag(fhom));
    hold on; plot(s.x,'.r');
    title('Exact Soln')
    plot(y_force.x,'.','MarkerSize',10,'LineWidth',10)
    quiver(real(y_force.x),imag(y_force.x),.2*pt_force(1:end/2),.2*pt_force(end/2+1:end))
    axis equal tight
end

N = 10*p*8;
[s,~] = quadr(s0, N, qtype, qntype); 
s = half_quadr(s);

% solve for tau using self close evaluation matrix
warning('off','MATLAB:nearlySingularMatrix')
if lptype == 's'    % Dirichlet
    f = AxiStokeslet(y_force,s)*pt_force;   % rhs
    A = AlpertSphereSLPMat(s);              % self eval matrix
    tau =  A\f;                             % density
else    % Neumann
    fp = AxiStokesletT(y_force,s)*pt_force;
    T = AlpertSphereSLPMatT(s);
    tau = (-eye(size(T))/2 + T)\fp;
end

fpres = AxiStokesletP(y_force,s)*pt_force;
Pres = SphereSLPMatP(s)*tau;

% evaluate velocity field
u = nan*(1+1i)*zz;
temp = SphereQuadtp(s,t)*tau;
u(ii) = temp(1:end/2) + 1i*temp(end/2+1:end);

 % evaluate pressure field
pres = nan*zz;
temp = SphereQuadtpP(s,t)*tau;
pres(ii) = temp;

% 
Nn = [Nn;numel(s.x)];
err = u(2*end/2:end,2*end/2:end)-fhom(2*end/2:end,2*end/2:end);
Err = [Err;max(abs(err(:)))];

% plot
figure(1),clf,imagesc(gx,gy,log10(abs(u-fhom))), axis equal
colorbar, hold on, plot(s.xlo,'or'), hold off, caxis([-12 0])    
title('log10 err in |u|')

figure(2),clf,imagesc(gx,gy,log10(abs(pres-fhompres))), axis equal
colorbar, hold on, plot(s.xlo,'or'), hold off, caxis([-12 0])    
title('log10 err in pressure')


keyboard

function G = AxiStokeslet(s,t)

X = [real(t.x);imag(t.x)]; Y = [real(s.x);imag(s.x)];
M = length(X)/2; N = length(Y)/2;

x1 = X(1:M); x2 = X(1+M:2*M);
y1 = Y(1:N); y2 = Y(1+N:2*N);  
G = zeros(2*M,2*N);  
for j =1:M
    Ker = AxiKernel(y1, y2, x1(j), x2(j));   
    G(j, 1:N) = ((Ker(:,2) + Ker(:,3)))';  
    G(j, N+1:end) = (Ker(:,4))';

    G(j+M, 1:N) = (Ker(:,5))';
    G(j+M, N+1:end) = ((Ker(:,1) + Ker(:,6)))';  
%     pause
end
W = 0.5*ones(size(s.x))/8/pi;  
G = G.*[repmat(W', 2*M, 1), repmat(W', 2*M, 1)];
G = 2*G;

end

function G = AxiStokesletT(s,t)

X = [real(s.x);imag(s.x)]; m = length(X)/2;
x1 = X(1:m); x2 = X(1+m:2*m);
% w = s.w;

mt = numel(t.x);
G1 = zeros(mt,m); G2 = G1; G3 = G1; G4 = G1;

for ne = 1:mt    % target
    
%     wt = w;
    wt = ones(size(s.x));
    Ker = AxiKernelT(x1, x2, real(t.x(ne)), imag(t.x(ne)), real(t.nx(ne)), imag(t.nx(ne))); 
    
    G1(ne, :) = (wt.*(Ker(:,1)))'; 
    G2(ne, :) = (wt.*Ker(:,2))';
    
    G3(ne, :) = (wt.*Ker(:,3))';
    G4(ne, :) = (wt.*(Ker(:,4)))';
end

W = -3/4/pi; 

G = [G1 G2;G3 G4]; 
G = G*W;


end

function G = AxiStokesletP(s,t)

X = [real(t.x);imag(t.x)]; Y = [real(s.x);imag(s.x)];
M = length(X)/2; N = length(Y)/2;

x1 = X(1:M); x2 = X(1+M:2*M);
y1 = Y(1:N); y2 = Y(1+N:2*N);  
G = zeros(M,2*N);  
for j =1:M
    Ker = AxiKernelP(y1, y2, x1(j), x2(j));   
    G(j, 1:N) = (Ker(:,1))';  
    G(j, N+1:end) = (Ker(:,2))';

end
W = ones(size(s.x))/4/pi;  
G = G.*[repmat(W', M, 1), repmat(W', M, 1)];
end
