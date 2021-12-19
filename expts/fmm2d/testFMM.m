% test FMM
% may need to copy libgfortran.so.3 and libgfortran.so.3.0.0 to /usr/local/MATLAB/R2021a/sys/os/glnxa64
% for recent release of Matlab...
% 
%
% Hai 01/20/21

addpath ..

clear s
a = .1; b = 8; 
s.Z = @(t) (1 + a*cos(b*t)).*exp(1i*t); s.p = 16; Np = 8*4;
[s, N, np] = quadr(s, s.p*Np, 'p', 'G');

% check SLP self
t.x = []; sig = s.cur;
[u1,~]=LapSLPFMM2D(s,1,t,0,sig);
A1 = LapSLPmat(s,s); A1(1:numel(s.x)+1:end) = 0;
u_dir1 = A1*sig(:);
diff1 = max(abs(u1(:)-u_dir1(:)))

% check SLP target
ntarget = 2000;
t.x = rand(1,ntarget)+1i*rand(1,ntarget); t.x = t.x(:);
[~,u2]=LapSLPFMM2D(s,0,t,1,sig);
A2 = LapSLPmat(t,s); 
u_dir2 = A2*sig(:);
diff2 = max(abs(u2(:)-u_dir2(:)))

% check DLP self
sig = ones(size(s.cur)); s.ws = ones(size(s.ws));
[u3,~]=LapDLPFMM2D(s,1,t,0,sig);
A3 = LapDLPmat(s,s); A3(1:numel(s.x)+1:end) = 0;
u_dir3 = A3*sig(:);
diff3 = max(abs(u3(:)-u_dir3(:)))

% check DLP target
[~,u4]=LapSLPFMM2D(s,0,t,1,sig);
A4 = LapSLPmat(t,s); 
u_dir4 = A4*sig(:);
diff4 = max(abs(u4(:)-u_dir4(:)))

keyboard


