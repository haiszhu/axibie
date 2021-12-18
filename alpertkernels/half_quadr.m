function s = half_quadr(s0)

s.Z=s0.Z; 
if isfield(s0,'Zp'), s.Zp = s0.Zp; end
if isfield(s0,'Zpp'), s.Zpp = s0.Zpp; end
s.p=s0.p; s.tlo=s0.tlo(1:end/2); s.thi=s0.thi(1:end/2); s.np=s0.np/2;
s.xlo=s0.xlo(1:end/2); s.xhi=s0.xhi(1:end/2); s.w=s0.w(1:end/2);
s.x=s0.x(1:end/2); s.xp=s0.xp(1:end/2); s.xpp=s0.xpp(1:end/2);
s.sp=s0.sp(1:end/2); s.tang=s0.tang(1:end/2); s.nx=s0.nx(1:end/2);
s.cur=s0.cur(1:end/2); s.w=s0.w(1:end/2); s.ws=s0.ws(1:end/2); s.t=s0.t(1:end/2); 
s.wxp=s0.wxp(1:end/2);

end