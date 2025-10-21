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
  persistent polar
  if nargin > 0 && ischar(varargin{1}) % set generator or query
    switch varargin{1}
      case 'Polar', polar = true;
      case 'BuiltIn', polar = false;
      case 'Query'
        if isempty(~polar) || ~polar
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
    if n<=2
      while true
        uv = rand_unif(2,1);
        u = uv(1)*2 - 1;
        v = uv(2)*2 - 1;
        s = u*u + v*v;
        if s > 0 && s < 1, break, end
      end
      R = sqrt(-2*log(s)/s);
      x(1) = u*R;
      if n==2, x(2) = v*R;
      end
    else
      k = 1;
      while k <= n
        uv = rand_unif(2,1);
        u = uv(1)*2 - 1;
        v = uv(2)*2 - 1;
        s = u*u + v*v;
        if s > 0 && s < 1
          R = sqrt(-2*log(s)/s);
          x(k) = u*R;
          k = k+1;
          if k <= n, x(k) = v*R; k = k+1; end
        end
      end
    end
    x = reshape(x, varargin{:});
  end
end
