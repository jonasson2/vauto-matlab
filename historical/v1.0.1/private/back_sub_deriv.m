% BACK_SUB_DERIV  Derivative of back substitution
%
%   Xd = BACK_SUB_DERIV(L, Ld, X, Yd) finds the derivative of the solution X to
%   an upper triangular system L'·X = Y. Ld and Yd should have derivatives of L
%   and Y, and Xd returns the derivative of X.
%
%   METHOD: From L'·X = Y it follows that Ld(:,:,i)'·X + L'·Xd(:,:,i) = Yd, and
%   thus Xd may be obtained with back substitution. 

function Xd = back_sub_deriv(L, Ld, X, Yd)
  LTT.LT = true; LTT.TRANSA = true;
  [n, m] = size(X);
  nPar = size(Ld,3);
  LdX = permute(reshape(reshape(Ld, [n,n*nPar])'*X, [n,nPar,m]), [1,3,2]);
  RHS = reshape(Yd - LdX, [n,m*nPar]);
  Xd = reshape(linsolve(L,RHS,LTT), [n,m,nPar]);
end
