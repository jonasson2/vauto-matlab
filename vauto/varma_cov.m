function C = varma_cov(X, maxlag, norm)

  % VARMA_COV  Sample autocovariance of multivariate time series X
  %
  %   C = VARMA_COV(X) computes the autocovariance matrices of a multivariate
  %   time series for all lags, with the t-th state, x(t), stored in the t-th
  %   column of X, t=1,...,n. C will be r by r by n with C(:,:,k) = Ck =
  %   autocovariance at lag k.
  %
  %   VARMA_COV(X, MAXLAG) computes C0,...,C{maxlag}.
  %
  %   VARMA_COV(X, MAXLAG, NORM) sets normalization, if NORM is "ML" (the
  %   default) Ck is normalized by 1/n, giving the maximum likelihood estimate
  %   of the theoretical autocovariance, and if NORM is "corr" or "corrected" Ck
  %   is normalized by 1/(n-k), giving an unbiased estimate for the case that
  %   the mean is known (e.g. known to be zero). NORM is not case sensitive.
  %
  %   Mathematically,
  %     Ck(i,j) = Cov(xi_t, xj_{t−k})
  %     Ck      = Cov(x_t, x_{t−k}) = Z·Y' * N,
  %
  %   where Y has the leading n−k columns of X, Z has the trailing n−k columns,
  %   and N is the normalization factor.

  [r, n] = size(X);
  if nargin < 3, norm = "ML"; end
  if nargin < 2, maxlag = size(X, 2) - 1; end
  if maxlag < 0 || maxlag >= n
    error('maxlag must satisfy 0 <= maxlag < n')
  end
  norm = lower(norm);
  if ~ismember(norm, ["ml", "corr", "corrected"])
    error("norm should be ML, corr or corrected");
  end
  mu = mean(X, 2);
  X = X - mu;
  C = zeros(r, r, maxlag + 1);
  for k = 0:maxlag
    if norm == "ml", f = n; else, f = n - k; end
    Y = X(:, 1+k:end);
    Z = X(:, 1:end-k);
    C(:, :, k+1) = Z*Y'/f;
  end
end
