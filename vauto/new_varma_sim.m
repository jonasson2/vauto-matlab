% NEW_VARMA_SIM  Simulate an ARMA or a VARMA time series model
% 
%   ARMA MODELS:
%     X = NEW_VARMA_SIM(A,B,sigma,n) with scalar sigma generates a zero-mean
%     time series of length n using the ARMA(p,q) model:
% 
%              x(t) = A1·x(t-1) + ... + Ap·x(t-p) + y(t)                    (1a)
%     where:
%              y(t) = eps(t) + B1·eps(t-1) + ... + Bq·eps(t-q)              (1b)
% 
%     and eps(t) is N(0,sigma). A should be the row vector [A1,...,Ap] and B
%     should be [B1,...,Bq]. X returns the series as a row vector. To start the
%     simulation, (x(1),...,x(h),eps(1),...,eps(h)), with h=max(p,q), are drawn
%     from N(0,SCE) where SCE = [SS CC; CC' EE] is obtained by solving the
%     modified Yule-Walker equations (see vyw_factorize, vyw_solve and S_build),
%     so there is no need to throw away the first values to avoid spin-up
%     effects. The model must be stationary. 
% 
%     X = NEW_VARMA_SIM(A,B,sigma,n,mu) returns a time series with mean mu using
%     y(t) as above, and x(t) given by:
% 
%             x(t) - mu = A1·(x(t-1) - mu) + ... + Ap·(x(t-p) - mu) + y(t)   (2)
% 
%     In this case x and eps are drawn from N([mu; 0], SCE).
% 
%     X = NEW_VARMA_SIM(A,B,sigma,n,mu,M) generates M sequences simultaneously
%     and returns the i-th one in the i-th column of X. Use mu=[] to obtain
%     zero-mean.
% 
%     X = NEW_VARMA_SIM(A,B,sigma,n,mu,M,X0) sets x(1)...x(h) to X0(end-h+1:end)
%     - mu and draws eps(1)...eps(h) from the conditional distribution of
%     eps(1), ..., eps(h) given x(1),...,x(h). For this option the model need
%     not be stationary.
% 
%     [X,EPS] = NEW_VARMA_SIM(...) returns also the shock series, eps(t), in the
%     t-th row of EPS.
% 
%   VARMA MODELS:
%     X = NEW_VARMA_SIM(A,B,Cov,n), where Cov is an r×r matrix uses a VARMA(p,q)
%     model given by (1), but now x(t) is r-dimensional, eps(t) is r-variate
%     normal with mean 0 and covariance Cov, and Cov and the Ai's and Bi's are
%     r×r matrices. A and B should contain A = [A1 A2...Ap] (an r × r·p matrix)
%     and B = [B1 B2 ... Bq] (an r × r·q matrix). The r×n matrix X returns x(t)
%     in its t-th column. As in the scalar case the simulated series is spin-up
%     free, the starting values x(1),...,x(h) and eps(1),...,eps(h) being drawn
%     from the correct distribution. The model must be stationary.
%  
%     X = NEW_VARMA_SIM(...,mu) uses (2) for x(t) instead of (1a).
% 
%     X = NEW_VARMA_SIM(...,mu,M) returns M sequences in an r×n×M
%     multidimensional X (when r > 1) with the i-th sequence in X(1:r, 1:n, i).
%     Use mu = [] to obtain zero-mean.
% 
%     X = NEW_VARMA_SIM(...,mu,M,X0) initializes x(1),...,x(h) with the last h
%     columns of X0 instead of drawing from N(mu,SS). The model need not be
%     stationary.
% 
%     [X,EPS] = NEW_VARMA_SIM(...) returns also the shock series; the i-th
%     eps(t) is returned in EPS(:, t, i).
% 
%   For both ARMA and VARMA, use NEW_VARMA_SIM(A,[],...) for a pure
%   autoregressive model, and NEW_VARMA_SIM([],B,...) for a pure moving average
%   model.
% 
%   The method used is described in [3]. It is an improved version of varma_sim
%   in the original Vauto package aka ACM TOMS Algorithm 878 as described in [1]
%   and [2]. The primary difference is that in this new function the first h
%   shocks (eps_t) are drawn first and afterwards the first h states (x_t) are
%   drawn using the conditional distribution of x_t|eps_t, whereas the original
%   version did the opposite (eps_t after x_t).
%
%   [1] K Jonasson and SE Ferrando 2008. Evaluating exact VARMA likelihood and
%       its gradient when data are incomplete. ACM Trans. Math. Softw. (TOMS),
%       35(1)
%
%   [2] K Jonasson 2008. Algorithm 878: Exact VARMA likelihood and its gradient
%       for complete and incomplete data with Matlab. ACM Trans. Math. Softw.
%       (TOMS), 35(1)
%
%   [3] K Jonasson 2025. Burn-in free simulation of VARMA time series.
%       Manuscript in preparation.
%
%   (C) Kristján Jónasson, Dept. of Computer Science, University of Iceland,
%   2025. jonasson@hi.is.

function [x, eps] = new_varma_sim(A, B, Sig, n, mu, M, x0)
  r = size(Sig, 1);
  if isempty(A), A = zeros(r,0); end
  if isempty(B), B = zeros(r,0); end
  [r, p, q, h] = getdims(A, B, Sig);
  Aflp = flipmat(A);
  Bflp = flipmat(B);
  if n<h, error('Too short series'); end
  if nargin < 5 || isempty(mu), mu = zeros(r,1); else mu = mu(:); end
  if nargin < 6 || isempty(M), M=1; end
  if nargin < 7, x0 = []; end
  x = zeros(n*r,M);
  eps = zeros(n*r,M);
  Psi = find_Psi(A, B);
  G = find_G(A, Bflp, Psi, Sig);
  PLU = vyw_factorize(A);
  assert(isempty(PLU) || isempty(PLU{1}) || PLU{1}(1) ~= 0)  % vyw_factorize ok
  S = vyw_solve(A, PLU, G);
  X = zeros(r*n, M);
  if isempty(x0)
    I = 1:r*h;
    SS = S_build(S, A, G, h); % SS = cov(x, x)
    Psi_hat = find_Psi_hat(Psi, Sig);
    R = SS - Psi_hat*Psi_hat';
    E = reshape(randnm(n*M, Sig)', [r*n, M]);
    X(I,:) = Psi*E(I,:) + randnm(M, R)';
  else  % x0 given
    % Find covariances and variances
    I = 1:r*k;
    k = length(x0);
    assert(h <= k && k <= n)
    CC = find_C(A, B, Sig, k);
    SS = S_build(S, A, G, k);

    % Find e = E(eps{1:k}|x0)
    X(I,:) = x0(:) - repmat(mu,k,1);
    E = zeros(r*n, M);
    e = CC'*LS'\(LS\X);

    % Find R = Var(eps{1:k}|x0)
    LS = chol(SS, 'lower');
    CC = LS\CC;
    R = -CC'*CC;
    J = 1:r;
    for j = 1:k
      R(J,J) = R(J,J) + Sig;
      J = J + r;
    end
    E(I,:) = e + randnm(M, R)';
    E(k*r+1:end) = reshape(randnm((n-k)*M, Sig)', [], M);
  end
  I = r*h + (1:r);
  J = r*(h-p)+1 : r*h-1;
  K = r*(h-q)+1 : r*h-1;
  for t = h+1:n
    X(I,:) = E(I,:) + Aflp*X(J,:) + Bflp*X(K,:);
    I = I+r;
    J = J+r;
    K = K+r;
  end
  if r==1 && M==1  %  one ARMA sequence:
    x = reshape(x,1,n) + mu;
    eps = reshape(eps,1,n);
  elseif r==1      %  several ARMA sequences:
    x = reshape(x,n,M) + mu;
    eps = reshape(eps,n,M);
  else             %   one or more VARMA sequences in r×n×M array:
    x = reshape(x,r,n,M) + repmat(mu,[1,n,M]);
    eps = reshape(eps,r,n,M); 
  end
end

function x = randnm(n, Sig)
  %  RANDNM  Multivariate normal random vectors
  r = size(Sig, 1);
  [R,p] = chol(Sig);
  assert(p==0);
  x = rand_norm(n, r)*R;
end

function Psi_hat = find_Psi_hat(Psi, Sig)
  LSig = chol(Sig, 'lower');
  r = size(Sig, 1);
  h = size(Psi, 1)/r;
  I = 1:r;
  Psi_hat = zeros(h, h);
  for k = 1:h
    Psi_hat(:,I) = Psi(:,I)*LSig;
    I = I+r;
  end
end
