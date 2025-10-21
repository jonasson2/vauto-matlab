%TEST_VARMA_LLC_DERIV
%
%  TEST_VARMA_LLC_DERIV compares gradient of log-likelihood function with a
%  numerical gradient for two testcases (the gradient is further checked with
%  test_omega_deriv).

function test_varma_llc_deriv(varargin)
  fprintf('TESTING VARMA_LLC_DERIV...')
  rand_init(1);
  quiet=0;
  if nargin==0, cases = num2cell(1:11); %{'mediumAR', 'mediumARMA1'};
  elseif ischar(varargin{1}), quiet = 1; cases = num2cell(1:11);
  else cases = num2cell(varargin{1}); end
  if ~quiet, fprintf('\n  Max analytic / numerical gradient difference:'); end
  for tcase = cases
    [A, B, Sig, p, q, r, name] = testcase(tcase{:});
    n = 10;
    X = varma_sim(A, B, Sig, n);
    theta = mat2theta(A, B, Sig);
    d = diff_test(@loglik, theta, X, p, q, r);
    if ~quiet, fprintf('  Testcase %-12s %.1e\n', [name ':'], d); end
    ascertain(d<1e-8)
  end
  disp('  OK')
end

function [ll,lld] = loglik(theta, X, p, q, r)
  [A, B, Sig] = theta2mat(theta, p, q, r);
  if nargout==1
    [ll, ok] = varma_llc(X, A, B, Sig);
  else
    [ll, ok, lld] = varma_llc(X, A, B, Sig);
  end
  ascertain(ok);
end
