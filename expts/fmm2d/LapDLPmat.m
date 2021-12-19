function A = LapDLPmat(t,s)
%
%
% Hai 01/20/21

d = bsxfun(@minus,t.x,s.x.');
A = real(bsxfun(@rdivide,(1/2/pi)*s.nx.',d));
A = bsxfun(@times, A, s.ws(:)');

end