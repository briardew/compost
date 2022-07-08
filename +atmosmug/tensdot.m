% TENSDOT  Compute the inner product of two tensors.
%
%    C = TENSDOT(A, B)  computes the inner product A*B of A and B by reshaping
%    them into matrices, computing their matrix product, and reshaping them
%    back.

% Author: Brad Weir
%==============================================================================%
function cc = tensdot(aa, bb)

nna = size(aa);
nnb = size(bb);

if (nna(end) ~= nnb(1))
  error('Size of last dimension of A must equal first dimension of B');
end

nn1 = prod(nna(1:end-1));
nn2 = nna(end);
nn3 = prod(nnb(2:end));

aam = reshape(aa,  nn1, nn2);
bbm = reshape(bb,  nn2, nn3);

ccm = aam * bbm;

cc  = reshape(ccm, [nna(1:end-1), nnb(2:end)]);
