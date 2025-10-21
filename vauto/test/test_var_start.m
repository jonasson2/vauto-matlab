%TEST_VAR_START  Test function var_start
%
%  TEST_VAR_START checks that the function var_start (which finds starting
%     values for likelihood maximization of VAR_LL) works correctly, by
%     comparing the output with that found with a more direct calculation). The
%     test is run for all but the last named testcase from "testcase.m" with
%     several different missing value patterns.
%  TEST_VAR_START(N) runs only the N-th named testcase.
%  TEST_VAR_START(N,M) runs N-th case on M-th missing value pattern.
%  TEST_VAR_START(...,'quiet') runs quietly.

function test_var_start(varargin)
  [args,quiet] = getflags(varargin,'quiet');
  fprintf('TESTING VAR_START');
  rand_init(1);
  if ~quiet, disp ' ', else fprintf('...'); end
  patt = [];
  if ~isempty(args)
    cases = {args(1)};
    if length(args) > 1, patt = args{2}; end
  else
    cases = num2cell(num2cell(1:testcase('number')-1));
  end
  ncases = length(cases);
  tol = 2e-14;
  for jcase = 1:ncases
    tcase = cases{jcase};
    [Ag, Bg, Sigg, p, q, r, name] = testcase(tcase{:});
    if q>0, continue; end
    mu = 0.1*(1:r)';
    tcstring = [int2str(tcase{1}) ' (' name ')'];
    fprintf_if(~quiet, 'TESTCASE %s\n', tcstring);
    r = size(Sigg,1);
    n = p + 10;
    Xfull = varma_sim(Ag, [], Sigg, n, mu);
    if ~isempty(patt), misspat = patt;
    else               misspat = [0 -4 -3 -2 -1 20];
    end
    for i=misspat
      casemsg = sprintf('testcase %d, missing pattern %d', jcase, i);
      fprintf_if(~quiet, '  missing pattern %d\n', i);
      X = makemissing(Xfull,i);
      miss = isnan(X);
      
      % Call function to be tested
      [A, Sig, mu] = var_start(X, p);

      % Fill missing values with "nanmeans" and subtract means from X-rows
      mu1 = zeros(r,1);
      for j = 1:r, mu1(j) = mean(X(j, ~miss(j,:))); end
      X = X - repmat(mu1, 1, n);
      X(isnan(X)) = 0;

      % Find conditional maximum likelihoods with "direct" linear regression
      Xhat = [];
      for j = 1:p, Xhat = [Xhat X(:,j:n-p-1+j)']; end
      bhat = X(:, p+1:n)';
      C = Xhat'*Xhat;
      b = Xhat'*bhat;
      SymPosDef = struct('SYM',true,'POSDEF',true);
      try
        A1 = linsolve(C, b, SymPosDef);
        tolA = tol;
      catch
        A1 = linsolve(C + 1e-10*eye(size(C)), b, SymPosDef);
        tolA = 1e-4;
      end
      resid = reshape(Xhat*A1 - bhat, n-p, r);
      A1 = A1';
      A1 = reshape(A1, r, r, p);
      A1 = flipdim(A1, 3);
      A1 = reshape(A1, r, r*p);

      % Compare results
      assert_almostequal(mu, mu1, tol, casemsg);
      assert_almostequal(A, A1, tolA, casemsg);
      assert_almostequal(tril(Sig), tril(cov(resid,1)), tol, casemsg);
    end
  end
  disp('  OK');
end
