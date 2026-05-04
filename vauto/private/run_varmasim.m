function run_varmasim(tc, n, M)
  if nargin < 3, M = 1; end
  if nargin < 2, n = 20; end
  rand_init('ParkMillerPolar', 42);
  [A, B, Sig] = testcase(tc);
  [X, ~] = new_varma_sim(A, B, Sig, n, 0, M);
  disp(X')
  mu = mean(X, 2);
  
end