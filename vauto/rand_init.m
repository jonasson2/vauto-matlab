
% RAND_INIT  Initialize random number generators
%
%   RAND_INIT('BuiltIn') sets both rand_unif and rand_norm to use the built-in
%   random number generators.
%
%   RAND_INIT('ParkMillerPolar') sets rand_unif to use the Park-Miller algorithm
%   and rand_norm to use the polar method. Allows getting same sequence as
%   programs written in C.
%
%   RAND_INIT(..., SEED) also sets the seed.

function rand_init(generator, seed)
  if isnumeric(generator)
    seed = generator;
    generator = 'BuiltIn';
  end
  switch generator
    case 'BuiltIn'
      rand_norm('BuiltIn');
      rand_unif('BuiltIn');
      rng('twister');
      if nargin > 1, rng(seed); end
    case 'ParkMillerPolar'
      rand_norm('Polar');
      rand_unif('ParkMiller');
      if nargin > 1
        rand_seed(seed); 
        rand_unif(); % spin_up one number
      end
    otherwise
      error('Illegal random number generator')
  end
end
