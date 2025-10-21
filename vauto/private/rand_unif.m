% RAND_UNIF  Uniform random numbers in [0,1)
%
%   Normally RAND_UNIF calls the built-in rand function of Matlab, but if the
%   call RAND_UNIF('ParkMiller') is made, it uses instead the Park-Miller linear
%   congruential generator described in [1]. This is provided so that one can
%   let programs written in C or Fortran produce exactly the same sequence of
%   random numbers as Matlab programs (to aid in testing and/or debugging). The
%   built-in rand function is about 1000 times faster.
%
%   X = RAND_UNIF() returns a uniformly distributed random number in [0,1)
%   X = RAND_UNIF(N) returns an N by N matrix of uniform random numbers.
%   X = RAND_UNIF(M,N) returns an M by N matrix of uniform random numbers.
%   X = RAND_UNIF(N1,N2,...) returns an N1 by N2 by ... array.
%
%   RAND_UNIF('ParkMiller') turns the ParkMiller generator on.
%   RAND_UNIF('BuiltIn') turns the built-in generator on.
%   RAND_UNIF('Query') returns the name of the active generator.
%
%   The built-in generator is set to the Mersenne Twister, being the default
%   generator in all versions since 2008b (until at least 2025b). It replaces
%   the SWB generator used in the original Vauto version in Vauto v1.1.0.
%
%   Use rand_init to set the state of the generator,
%
%   Reference
%   [1] S.K. Park and K.W. Miller (1988). "Randomn number generators: Good ones
%       are hard to find". Communications of the ACM 31 (10), 1192-1201.

function x = rand_unif(varargin)
  persistent a m q r fct
  if nargin > 0 && ischar(varargin{1}) % set generator or query
    switch varargin{1}
      case 'ParkMiller'
        if isempty(a)
          a = int32(16807);
          m = int32(2147483647);
          q = idivide(m,a);
          r = mod(m,a);
          fct = 1/(double(m)+1);
        end
      case 'BuiltIn'
        a = [];
      case 'Query'
        if isempty(a), x = 'BuiltIn';
        else, x = 'ParkMiller'; end
        return
      otherwise
        error('Illegal keyword')
    end
  else % generate random number
    if isempty(a) % built-in
      x = rand(varargin{:});
    else % Park-Miller
      seed = rand_seed();
      if nargin == 0, varargin = {1,1}; end
      if nargin == 1, varargin{2} = varargin{1}; end
      n = prod(cell2mat(varargin));
      x = zeros(n,1);
      for i=1:n
        s1 = a*mod(seed,q);
        s2 = r*idivide(seed,q);
        if s1 > s2, seed = s1 - s2;
        else        seed = (m - s2) + s1; end
        x(i) = double(seed)*fct;
      end
      x = reshape(x, varargin{:});
      rand_seed(seed);
    end
  end
end
