%

clear all
addpath(genpath('../'));
setup();

side = 'i'; % test interior or exterior
lptype = 's'; % test SLP or DLP
v = 1;
a = .3; w = 5; R = @(t) (1 + a*cos(w*t))*1; s.Z = @(t) R(t).*exp(1i*t); 
s.p = 16; N = s.p*64; 
% s.p = 16; N = s.p*40; % why only this case works for DLP... this is 
s = quadr(s, N, 'p', 'G');

% target
nx = 151; gx = ((1:nx)/nx*2-1)*1.5; ny = 151; gy = ((1:ny)/ny*2-1)*1.5; % set up plotting grid
[xx, yy] = meshgrid(gx,gy); zz = (xx+1i*yy);

% generate the exact solution
alpha = linspace(-pi/8,pi/8,5); y_force = [];
if side == 'e', y_force.x = 0.4*exp(1i*alpha).'; else, y_force.x = 2*exp(1i*alpha).'; end  % location of the point forces 
pt_force = [1;1;-1;1;-1/2];  % magnitude & direction of the point forces (sample 2)
y_force.ws = ones(size(y_force.x));

fhom = nan*xx; % exact soln
t = [];%s.x = [0;s.x];
[IN, ON] = inpolygon(real(zz),imag(zz),real(s.x),imag(s.x)); 
if side == 'e', ii = ~IN&~ON; else, ii = IN&~ON; end

t.x = zz(ii(:));  % eval pts only on one side
A = LapSLPmat(t,y_force);
f_temp = A*pt_force;
fhom(ii(:)) = f_temp; % the exact soln

% plot the exact soluntion and the point forces
if v==1
    figure(2), clf, streamslice(gx,gy,real(fhom),imag(fhom));
    contourf(xx,yy,fhom,10)
    hold on; plot(real(s.Z(s.t)),imag(s.Z(s.t)),'r');
    title('Exact Soln')
    plot(y_force.x,'.','MarkerSize',10,'LineWidth',10)
    axis equal tight
end

% bdry condition
f = LapSLPmat(s,y_force)*pt_force; 

% Laplace self interaction matrix
if lptype == 's'
    A = LapSLPSpecialMatDG(s,s);
end
if lptype == 'd'
    if side == 'i'
        A = -1/2*eye(numel(s.x))+LapDLPSpecialMatDG(s,s,1); % on-surface flag for DLP
    end
    if side == 'e'
        A = 1/2*eye(numel(s.x))+LapDLPSpecialMatDG(s,s,1) + LapSLPSpecialMatDG(s,s);
    end
%     A = LapDLPSpecialMat(s,s);
end

% solve
tau = A\f;

% eval at target
if lptype == 's'
    uLap = LapSLPSpecialMatDG(t,s)*tau;
end
if lptype == 'd'
    if side == 'i'
        uLap = LapDLPSpecialMatDG(t,s)*tau;
    end
    if side == 'e'
        uLap = (LapSLPSpecialMatDG(t,s)+LapDLPSpecialMatDG(t,s))*tau;
    end
end


% solution
u = nan*(1+1i)*zz;
u(ii(:)) = uLap;

% error
figure(2),clf,imagesc(gx,gy,log10(abs(u-fhom))),colorbar
hold on, axis equal, plot(s.x,'r-')
caxis([-15 -10]), colormap('jet'), title('error')

% check traction...

    
keyboard