%TEST_VAR_LL  Test the pure autoregressive likelihood function var_ll
%
%  TEST_VAR_LL compares likelihood calculated with var_ll with that calculated
%              with varma_llm for test cases that "testcase" creates. The
%              gradient is also compared, along with some further tests.
%  TEST_VAR_LL(N) runs only the N-th named testcase from "testcase"
%  TEST_VAR_LL(N,M) runs N-th case on M-th missing value pattern
%  TEST_VARMA_LLM FULL runs a thorough test with more combinations of (p,r) and
%                      more missing value patterns

function test_var_ll(varargin)
  [args,quiet] = getflags(varargin,'quiet');
  [args,full] = getflags(args,'full');
  fprintf('TESTING VAR_LL'); if full, fprintf(' (COMPREHENSIVE CHECK)'); end
  rand_init(1);
  if ~quiet, disp ' ', else fprintf('...'); end
  patt = [];
  if full
    cases = {};
    for p=0:4, for r=1:4, cases{end+1} = {p,0,r}; end, end
    tol = [5e-15 5e-14];
  elseif ~isempty(args)
    cases = {args(1)};
    if length(args) > 1, patt = args{2}; end
    tol = [5e-16 5e-15];
  else
    cases = num2cell(num2cell(1:testcase('number')-1));
    tol = [5e-16 5e-15];
  end
  % Check non-stationary models:......TODO
  A = [0.5 0.25];
  X = zeros(1,10); miss = isnan(X);
  mu = 0;
  msg = ['An indefinite Sig should return ok = 0 (false) and a positive '...
    'definite Sig should return ok = 1 (true)'];
  Sig = -1; [l, ok] = var_ll(X, A, Sig, mu, miss); assert_equal(ok, 0, msg);
  Sig = 1;  [l, ok] = var_ll(X, A, Sig, mu, miss); assert_equal(ok, 1, msg);
  ncases = length(cases);
  fnDiffMax = 0; gradDiffMax = 0;
  for j = 1:ncases
    tcase = cases{j};
    [A, B, Sig, p, q, r, name] = testcase(tcase{:});
    if full, tcstring = sprintf('p=%d q=%d r=%d', tcase{:});
    else tcstring = [int2str(tcase{1}) ' (' name ')']; end
    if q==0
      fprintf_if(~quiet, 'TESTCASE %s\n', tcstring);
      mu = 0.01*(1:r)';
      r = size(Sig,1);
      n = p + 6;
      X = varma_sim(A, [], Sig, n, mu);
      compare(A, Sig, X);
      if ~isempty(patt), misspat = patt;
      elseif ~full,      misspat = [0 -4 -3 -2 -1 20];
      elseif  full,      misspat = [0 -4 -3 -2 -1 2 4 10 15 20 25]; 
      end
      for i=misspat
        casemsg = sprintf('\nTest-case %d, missing pattern %d.', j, i);
        fprintf_if(~quiet, '  Missing pattern %d:', i);
        X1 = makemissing(X,i);
        miss = isnan(X1);
        %
        % COMPARE FUNCTION VALUES WITH THOSE OF VARMA_LLM
        [l1, ok] = varma_llm(X1, A, [], Sig, mu, miss);
        [l, ok] = var_ll(X1, A, Sig, mu, miss); assert_equal(ok, 1)
        msg = 'comparing function values from var_ll and varma_llm';
        fndiff = assert_equal(l, l1, 5e-14, [msg char(10) casemsg]);
        %
        % CHECK THAT LOGDETL IS OK WHEN RES IS ALSO RETURNED
        [l2, ok, res, xm] = var_ll(X1, A, Sig, mu, miss, 'res_miss');
        msg = 'Returning res and xm should not affect log-likelihood value';
        assert_equal(l, l2, [msg char(10) casemsg]);
        %
        % CHECK THAT RESIDUAL AND MISSING VALUE ESTIMATES ARE AS FOR VARMA_LLM
        [l2, ok, res1, xm1] = varma_llm(X1, A, [], Sig, mu, miss, 'res_miss'); 
        assert_equal(ok, 1);
        if ~isempty(xm1), assert_equal(xm, xm1, 1e-13); end
        %
        % COMPARE GRADIENTS
        [l, ok, ld] = var_ll(X1, A, Sig, mu, miss);
        [l, ok, ld1] = varma_llm(X1, A, [], Sig, mu, miss);
        msg =  ['ld has wrong length' char(10) casemsg];
        assert_equal(length(ld), p*r^2+r*(r+1)/2+r, msg);
        [n1, n2] = deal('varma_llm gradient', 'var_ll gradient');
        graddiff = assert_equal(ld, ld1, 5e-14);
        %
        % PRINT LINE TO REPORT
        fmt = ' max fn.diff.=%.1e, max grad.diff.=%.1e\n';
        fprintf_if(~quiet, fmt, fndiff, graddiff);
        %
        % UPDATE MAX DIFFERENCES
        fnDiffMax = max(fnDiffMax, fndiff);
        gradDiffMax = max(gradDiffMax, graddiff);
        %
        % CHECK THAT EMPTY MU RETURNS THE SAME AS ZERO MU:
        smallcase = p < 3 && q < 3 && r < 4;
        if smallcase
          zmu = zeros(r,1);
          [l0, ok0] = var_ll(X1, A, Sig, zmu, miss);
          [le, oke] = var_ll(X1, A, Sig, [], miss);
          msg=['Zero mu and empty mu should give same function value' casemsg];
          assert_equal(ok0, 1, msg);
          assert_equal(oke, 1, msg);
          assert_equal(l0, le, msg);
          [le, oke, lde] = var_ll(X1, A, Sig, [], miss);
          [l0, ok0, ld0] = var_ll(X1, A, Sig, zmu, miss);
          assert_equal(ok0, 1);
          assert_equal(oke, 1);
          assert_equal(length(lde), p*r^2+r*(r+1)/2, casemsg)
          assert_equal(length(ld0), p*r^2+r*(r+1)/2+r, casemsg)
          assert_equal(lde, ld0(1:p*r^2+r*(r+1)/2), 5e-14, casemsg)
          assert_equal(ok, 1, casemsg);
        end
      end;
    end
  end
  fmt = '  Overall, max function value diff.=%.1e, max gradient diff.=%.1e.\n';
  fprintf_if(~quiet, fmt, fnDiffMax, gradDiffMax);
  s1 = '  ("Max fn.diff." is difference between var_ll and varma_llm function';
  s2 = '  values and "max grad.diff." between gradients. The differences are';
  s3 = '  relative if max absolute value is > 1, otherwise they are absolute.)';
  fprintf_if(~quiet, '%s\n%s\n', s1, s2, s3)
  disp('  OK');
end

function compare(A, Sig, X)
  % Compare function values of (a) var_ll with zero mu and all-false miss, (b)
  % var_ll without mu- and miss-parameters and (c) varma_llc. Also check
  % redidual calculation.
  r = size(Sig,1);
  miss = false(size(X));
  ll = varma_llc(X, A, [], Sig);
  ll1 = var_ll(X, A, Sig, zeros(r,1), miss);
  ll2 = var_ll(X, A, Sig);
  [ll3, ok, res] = varma_llm(X, A, [], Sig, zeros(r,1), miss, 'res');
  [ll4, ok, res1] = var_ll(X, A, Sig, zeros(r,1), miss, 'res');
  [ll5, ok, res2] = var_ll(X, A, Sig, 'res');
  assert_equal(ll1, ll, 5e-14);
  assert_equal(ll2, ll, 5e-14);
  assert_equal(ll3, ll, 5e-14);
  assert_equal(ll4, ll, 5e-14);
  assert_equal(ll5, ll, 5e-14);
  assert_equal(res1, res, 5e-14);
  assert_equal(res2, res, 5e-14);
end
