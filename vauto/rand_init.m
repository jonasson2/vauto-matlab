% RAND_INIT  Initialize random number generators
%
%   RAND_INIT('BuiltIn') sets both rand_unif and rand_norm to use the built-in
%   random generators.
%
%   RAND_INIT('ParkMillerPolar') sets rand_unif to use the Park-Miller
%   algorithm and rand_norm to use the polar method. Used to be able to get same
%   sequence as programs written in C.
%
%   RAND_INIT(S) sets the state of both rand_unif and rand_norm to S.

function rand_init(s)
  if ischar(s)
    switch s
      case 'BuiltIn'
        rand_norm('BuiltIn');
        rand_unif('BuiltIn');
      case 'ParkMillerPolar'
        rand_norm('Polar');
        rand_unif('ParkMiller');
      otherwise
        error('Illegal keyword')
    end
  else
    rng(s, 'twister') % Changed from rand('state', s') and randn('state', s)
    if s < 100000, rand_seed(123456*s); else rand_seed(s); end;
  end
end
