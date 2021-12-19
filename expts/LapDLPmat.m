function A = LapDLPmat(t,s)
%
% modified based on the velocity kernel from BIE2D package, 12/26/20  

d = bsxfun(@minus,t.x,s.x.');                 % C-# displacements mat
% mult by identical rows given by source normals...
A = real(bsxfun(@rdivide,(1/2/pi)*s.nx.',d)); % complex form of dipole
A = bsxfun(@times, A, s.ws(:)');

end
