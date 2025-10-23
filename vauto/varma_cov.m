% VARMA_COV  Compute theoretical covariance of an ARMA or a VARMA time series
%
%   S0 = VARMA_COV(A, B, Cov) computes the theoretical variance of an
%   r-dimensional VARMA(p,q) model:
%
%     x(t) - mu = A1·(x(t-1) - mu) + ... + Ap·(x(t-p) - mu) + y(t)
%
%   where y(t) = eps(t) + B1·eps(t-1) + ... + Bq·eps(t-q) and eps is N(0,Cov).
%   A should be [A1...Ap] and B should be [B1...Bq]. As S0 is independent of mu
%   it is not needed as a parameter. When p=0 A can be empty and when q=0 B can
%   be empty.
%
%   S = VARMA_COV(A, B, Cov, k] computes a cell matrix with the covariances at
%   lags 0,1,...,k.
%
%   See varma_sim for references and description and discussion of the models.

function S = varma_cov(A, B, Sig, k)
  r = size(Sig, 1);
  if isempty(A), A = zeros(r,0); end
  if isempty(B), B = zeros(r,0); end
  p = size(A, 2)/r;
  if nargin < 4, k = 0; end
  [~, G] = find_CGW(A, B, Sig);
  PLU = vyw_factorize(A);
  if ~isempty(PLU) && ~isempty(PLU{1}) && PLU{1}(1) == 0
    error('Non-stationary model');
  end
  S = vyw_solve(A, PLU, G);
  if nargin < 4
    S = S{1};
  else
    A = reshape(A, r, r, p);
    for j = p+2:k+2
      S{j} = zeros(r,r);
      for i = 1:p
        S{j} = S{j} + A(:, :, i)*S{j-i};
      end
    end
  end
end
