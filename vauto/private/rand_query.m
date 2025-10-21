% RAND_QUERY  Return currently active random number generators
%
%   [UNIF, NORM] = RAND_QUERY returns the name of the active uniform RNG as
%   ParkMiller or BuiltIn and the active normal RNG as Polar or BuiltIn.

function [unif, norm] = rand_query()
  unif = rand_unif('Query');
  norm = rand_norm('Query');
end
