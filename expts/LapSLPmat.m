function A = LapSLPmat(t,s)
%
%
% modified based on the velocity kernel from BIE2D package, 12/26/20  

d = bsxfun(@minus,t.x,s.x.'); 
A = bsxfun(@times, log(abs(d)), -(1/2/pi)*s.ws(:)');  

end
