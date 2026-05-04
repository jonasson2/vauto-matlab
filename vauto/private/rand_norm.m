% RAND_NORM  Standard normal random numbers
%
%   Normally RAND_NORM calls the built-in randn function of Matlab, but if 
%   RAND_NORM('Polar') is called, it uses instead the polar method and with
%   input from RAND_UNIF with the Park-Miller method. This is provided for so
%   that one can let programs written in C or Fortran produce exactly the same
%   sequence of random numbers as Matlab programs (to aid in testing and/or
%   debugging). The built-in randn function is much faster.
%
%   X = RAND_NORM() returns an N(0,1) random number
%   X = RAND_NORM(N) returns an N by N matrix of N(0,1) random numbers.
%   X = RAND_NORM(M,N) returns an M by N matrix of N(0,1) random numbers.
%   X = RAND_NORM(N1,N2...) returns an N1 by N2 by ... multidimensional array.
%   RAND_NORM('Polar') lets future RAND_NORM calls use polar method.
%   RAND_NORM('BuiltIn') swithches back to built-in randn.
%   RAND_NORM('Query') returns the name of the active generator.
%
%   Use RAND_INIT to set the state for RAND_NORM.

function x = rand_norm(varargin)
  persistent polar spare_norm
  if nargin > 0 && ischar(varargin{1}) % set generator or query
    switch varargin{1}
      case 'Polar', polar = true; spare_norm = NaN;
      case 'BuiltIn', polar = false; spare_norm = NaN;
      case 'Query'
        if isempty(polar) || ~polar
          x = 'BuiltIn';
        else
          x = 'Polar';
        end
        return
      otherwise, error('Illegal keyword');
    end
  elseif isempty(polar) || ~polar % built-in
    x = randn(varargin{:});
else % polar
  if nargin == 0, varargin = {1,1}; end
  if nargin == 1, varargin{2} = varargin{1}; end
  n = prod([varargin{:}]);
  x = zeros(n,1);

  % initialise spare_norm the first time
  if isempty(spare_norm)
    spare_norm = NaN;
  end

  k = 1;

  % use cached spare if available
  if ~isnan(spare_norm)
    x(k) = spare_norm;
    spare_norm = NaN;
    k = k+1;
  end

  % fill as many pairs as possible
  while k+1 <= n
    while true
      uv = rand_unif(2,1);
      u = uv(1)*2 - 1;
      v = uv(2)*2 - 1;
      s = u*u + v*v;
      if s > 0 && s < 1
        break
      end
    end
    R = sqrt(-2*log(s)/s);
    x(k)   = u*R;
    x(k+1) = v*R;
    k = k+2;
  end

  % if one value left, generate a pair and cache the spare
  if k <= n
    while true
      uv = rand_unif(2,1);
      u = uv(1)*2 - 1;
      v = uv(2)*2 - 1;
      s = u*u + v*v;
      if s > 0 && s < 1
        break
      end
    end
    R = sqrt(-2*log(s)/s);
    x(k) = u*R;
    spare_norm = v*R;
  end

  x = reshape(x, varargin{:});  end
end
