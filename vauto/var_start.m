%   FIND STARTING VALUES FOR VAR LIKELIHOOD MAXIMIZATION.
%   
%   [A,B,Sig,mu] = VAR_START(X, p) determines A = {A1 A2...Ap}, Sig and mu that
%   can be used as starting values for numerical maximization of a VAR
%   likelihood function. The r×n array X contains the observed time series with
%   NaN in missing value positions; p is the number of autoregressive terms.
%   
%   METHOD: The Ai-s are chosen to minimize the residual sum of squares (or
%   equivalently maximize the conditional likelihood). Sig is chosen as the data
%   covariance matrix of the residuals. If there are missing values, these are
%   first filled in with the average of the corresponding series (row in X).

function [A, Sig, mu] = var_start(X, p)
  [r,n] = size(X);
  mu = zeros(r,1);
  miss = isnan(X);
  anyobs = any(~miss,2);
  for i = 1:r
    if anyobs(i), mu(i) = mean(X(i, ~miss(i,:))); end
    X(i,miss(i,:)) = mu(i);   % fill in missing values
  end
  X = X - repmat(mu, 1, n);
  [r, n] = size(X);
  x = X(:);
  if p > 0
    N = p*r^2;
    xd = zeros(r*n, 1, N);
    A = zeros(r,r*p);
    F = zeros(r,r,p,r,r,p);
    b = zeros(r,r,p);
    G = zeros(r,r,p+1);
    for d=0:p
      G(:,:,d+1) = X(:,p+1-d:n-p)*X(:,p+1:n-p+d)';
    end
    for i=0:p
      for j=i:p
        d = j-i; ne = n-p+1;
        V = G(:,:,d+1) + X(:,p+1-j:p-d)*X(:,p+1-i:p)'+X(:,ne:n-j)*X(:,ne+d:n-i)';
        for l=1:r
          if i>0,  F(l,:,j,l,:,i) = V; end
          if i==0 && j>0, b(l,:,j) = V(:,l); end
        end
      end
    end
    F = reshape(F,N,N); b = reshape(b,N,1);
    SymPosDef = struct('SYM',true,'POSDEF',true);
    ok = false; del = 1e-10;
    warning off MATLAB:nearlySingularMatrix
    while ~ok
      try
        a = linsolve(F', b, SymPosDef);
        ok = true;
      catch
        F = F + del*eye(N);
        del = del*10;
      end
    end
    warning on MATLAB:nearlySingularMatrix
    A = reshape(a,r,r*p);
    w = lambda_multiply(A, X(:), false(r, n));
    W = reshape(w(r*p+1:end), r, n-p);
    Sig = cov(W', 1); % normalize with n-p
  else  % No autoregressive terms. Find pairwise covariances
    A = zeros(r, 0);
    Sig = zeros(r);
    for i=1:r
      for j=1:i
        pij = ~miss(i,:) & ~miss(j,:);
        nij = sum(pij);
        Sig(i,j) = X(i,pij)*X(j,pij)'/nij;
      end
    end
    Sig = Sig + tril(Sig,-1)';
    [R,pidx] = chol(Sig);
    if pidx~=0, %  Add delta to diagonal of Sig
      del = max(1,norm(diag(Sig),inf))*0.0001;
      while pidx~=0
        Sig = Sig + del*eye(r);
        [R,pidx] = chol(Sig);
        del = del*10;
      end
    end
  end
end
