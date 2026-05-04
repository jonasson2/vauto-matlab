% IS_STATIONARY  Check VARMA model for stationarity
%   IS_STATIONARY(A) returns true if the VARMA-process x(t) = A1·x(t-1) + ... +
%   Ap·x(t-p) + y(t) is stationary, false otherwise. A is the r × r·p matrix
%   [A1...Ap]. The time series y(t) is given by a moving average process, y(t) =
%   B1·y(t-1) + ... + Bq·y(t-q) + eps(t), where eps(t) is N(0,Sig). Whether x is
%   stationary or not does not depend on B, so it need not be a parameter.
%
%   IS_STATIONARY(A, PLU) uses PLU from a previous call to vyw_factorize.
%
%   [STAT, Su] = IS_STATIONARY(...) returns also in Su the covariance matrix of
%   [x(t);...;x(t+p-1)], which is positive definite iff the process is
%   stationary.
%
%   IS_STATIONARY(A, 'specrad') uses the spectral radius of the companion
%   matrix, rho(Ac) where:
%
%       Ac = A1  A2  ... Ap-1  Ap
%             I
%                 I
%                    ...
%                          I    O
%
%   which measures the asymptotic growth rate of the process excluding shocks,
%   the process being stationary iff rho < 1.
%
%   [STAT, RHO] = IS_STATIONARY(..., 'specrad') returns rho.
%
%   In the non-specrad case, the test is based on solving the modified
%   Yule-Walker equations for the AR process x(t)=A1·x(t-1)+...+Ap·x(t-p)+eps(t)
%   with eps(t) ~ N(0,I), but the solution to these equations gives a postive
%   definite Su matrix (see omega_build) if and only if all the roots of the
%   polynomial fi(b) = I - A1·b - A2·b^2 - ... - Ap·B^p are outside the unit
%   circle, which characterizes a stationary VARMA process.

function [stat, argout] = is_stationary(A, argin)
  if nargin == 2 && (ischar(argin) || isstring(argin)) && strcmp(argin, 'specrad')
    if isempty(A)
      rho = 0;
    else
      r = size(A, 1);
      p = length(A)/r;
      k = r*(p-1);
      Ac = [A; eye(k) zeros(k,r)];
      rho = max(abs(eig(Ac)));
    end
    stat = rho < 1;
    argout = rho;
  else
    if isempty(A), stat=true; return, end
    r = size(A,1);
    p = length(A)/r;
    if nargin==1, PLU = vyw_factorize(A); end
    if nargin==2, PLU = argin; end
    if ~isempty(PLU) && ~isempty(PLU{1}) && PLU{1}(1) == 0
      stat=false;
      return
    end
    I = eye(r);
    S = vyw_solve(A, PLU, {I});
    Su = omega_build(S, {I}, {I}, p, p);
    Su = tril(Su) + tril(Su,-1)';
    [~, pp] = chol(Su);
    stat = pp==0;
    argout = Su;
  end
end
