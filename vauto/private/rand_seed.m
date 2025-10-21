% RAND_SEED  Utility to store random number seed
%
%   RAND_SEED(S) stores the seed S in a persistent variable.
%   S = RAND_SEED retrieves the seed previously stored.
%
% Used by the gateway functions in rand_norm.c and rand_unif.c

function sout = rand_seed(s)
  persistent seed
  if isempty(seed), seed = 123456; end
  if nargin > 0, seed = s; end
  if nargout > 0, sout = seed; end
end
