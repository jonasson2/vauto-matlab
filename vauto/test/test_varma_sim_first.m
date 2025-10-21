% TEST_VARMA_SIM_FIRST  First test of varma_sim.
%
%   TEST_VARMA_SIM calls varma_sim for a few simple cases and does a few
%   rudimentary checks on the results.
%
%   A (much) more comprehensive test of varma_sim is done by test_varma_sim

function test_varma_sim_first
  fprintf('TESTING VARMA_SIM (PRELIMINARY TEST)... ');
  for i=1:6
    [A, B, Sig, p, q, r, name] = testcase(i);
    h = max(p,q);
    %
    X1 = varma_sim(A,B,Sig,10);
    X2 = varma_sim(A,B,Sig,10,[],1);
    X3 = varma_sim(A,B,Sig,10,[],5); 
    [X4,eps4] = varma_sim(A,B,Sig,10);
    [X5,eps5] = varma_sim(A,B,Sig,10,[],5); 
    assertsize(X1,r,10);
    assertsize(X2,r,10);
    assertsize(X4,r,10); assertsize(eps4,r,10); 
    if r==1,
      assertsize(X3,10,5);
      assertsize(X5,10,5); assertsize(eps5,10,5); 
    else
      assertsize(X3,r,10,5);
      assertsize(X5,r,10,5); assertsize(eps5,r,10,5); 
    end
  end
  fprintf('OK\n');
end

function assertsize(A,n,m,k)
  if nargin == 3 || k==1
    ascertain(isequal(size(A), [n,m]));
  else
    ascertain(isequal(size(A), [n,m,k]));
  end
end
