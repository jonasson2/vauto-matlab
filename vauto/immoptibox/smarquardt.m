function  [X,info,perf,B] = smarquardt(fun, x0, opts, B0, varargin)
%SMARQUARDT  Secant version of Levenberg-Marquardt's method for least 
% squares: Find  xm = argmin{F(x)} , where  x = [x_1, ..., x_n]  and
% F(x) = .5 * sum(f_i(x)^2) .
% The functions  f_i(x) (i=1,...,m)  must be given by a MATLAB
% function with declaration
%            function  f = fun(x, p1,p2,...)
% p1,p2,... are parameters of the function.  In connection with 
% nonlinear data fitting they may be arrays with coordinates of 
% the data points.
%
% Call
%   [X, info] = smarquardt(fun, x0)
%   [X, info] = smarquardt(fun, x0, opts)
%   [X, info] = smarquardt(fun, x0, opts, B0, p1,p2,...)
%   [X, info, perf] = smarquardt(.....)
%   [X, info, perf, B] = smarquardt(.....)
%
% Input parameters
% fun  :  Handle to the function.
% x0   :  Starting guess for  xm .
% opts :  Vector with at most five elements.  
%         opts(1)  used in starting value for Marquardt parameter: 
%             mu = opts(1) * max(A0(i,i))  with  A0 = J(x0)'*J(x0)
%         opts(2:4)  used in stopping criteria:
%             ||F'||inf <= opts(2)                     or 
%             ||dx||2 <= opts(3)*(opts(3) + ||x||2)    or
%             no. of iteration steps exceeds  opts(4) .
%         opts(5)  "relative" step length for difference approximations.
%         Default  opts = [1e-3 1e-4 1e-8 100 1e-7]
%         If the input opts has less than 5 elements, it is
%         augmented by the default values.  
% B0   :  Initial approximation to J.
%         If  B0 is not given, a forward difference approximation
%         to it is used.
% p1,p2,..  are passed dirctly to the function FUN .    
%
% Output parameters
% X    :  If  perf  is present, then array, holding the iterates
%         columnwise.  Otherwise, computed solution vector.
% info :  Performance information, vector with 7 elements:
%         info(1:4) = final values of 
%             [F(x)  ||F'||inf  ||dx||2  mu/max(A(i,i))] ,
%           where  A = B'* B .
%         info(5) = no. of iteration steps
%         info(6) = 1 :  Stopped by small gradient
%                   2 :  Stopped by small x-step
%                   3 :  No. of iteration steps exceeded 
%                  -1 :  x is not a real valued vector
%                  -2 :  f is not a real valued column vector
%                  -4 :  Dimension mismatch in x, f, B0
%                  -5 :  Overflow during computation
%                  -6 :  Error in approximate Jacobian
%         info(7) = no. of function evaluations
% perf :  Array, holding 
%         perf(1,:) = values of  F(x)
%         perf(2,:) = values of  || B'*f(x) ||inf
%         perf(3,:) = mu-values.
% B    :  Computed approximation to Jacobian at the solution.
%
% Method
% Approximate Gauss-Newton with Levenberg-Marquardt damping and 
% successive updating of Jacobian approximation.

% Version 04.04.10.  hbn(a)imm.dtu.dk

% Check parameters and function call
F = NaN;  ng = NaN;
info = zeros(1,7);
if  nargin < 2,  stop = -1;
else
  [stop x n] = checkx(x0);   
  if  ~stop
    [stop F f] = checkfJ(fun,x0,varargin{:});  info(7) = 1;
    if  ~stop
      %  Finish initialization
      if  nargin < 3,  opts = []; end
      opts = checkopts(opts, [1e-3 1e-4 1e-8 100 1e-7]); 
      % Jacobian
      if  nargin > 3  % B0 is given
        sB = size(B0);  m = length(f);
        if  sum(sB) == 0  % placeholder
          [stop B] = Dapprox(fun,x,opts(5),f,varargin{:});  
          info(7) = info(7) + n;
        elseif  any(sB ~= [m n]),  stop = -4;
        else,                      B = B0;   end
      else
        [stop B] = Dapprox(fun,x,opts(5),f,varargin{:});  
        info(7) = info(7) + n;
      end
      % Check gradient and J'*J
      if  ~stop
        g = B'*f;   ng = norm(g,inf);  A = B'*B;
        if  isinf(ng) | isinf(norm(A(:),inf)),  stop = -5; end 
      end
    end
  end
end
if  stop
  X = x0;  perf = [];  info(1:6) = [F ng  0  opts(1)  0  stop];
  return
end

% Finish initialization
mu = opts(1) * max(diag(A));    kmax = opts(4);
Trace = nargout > 2;
if  Trace
  X = repmat(x,1,kmax+1);
  perf = repmat([F; ng; mu],1,kmax+1);
end 


% Iterate
k = 1;   nu = 2;   nh = 0;
ng0 = ng;
ku = 0;  % direction of last update

while  ~stop
  if  ng <= opts(2),  stop = 1; 
  else
    [h mu] = geth(A,g,mu);
    nh = norm(h);   nx = opts(3) + norm(x);
    if  nh <= opts(3)*nx,  stop = 2; end 
  end 
  if  ~stop
    xnew = x + h;    h = xnew - x;  
    [stop Fn fn] = checkfJ(fun,xnew,varargin{:});  info(7) = info(7)+1;
    if  ~stop
      % Update  B
      ku = mod(ku,n) + 1; 
      if  abs(h(ku)) < .8*norm(h)  % extra step
        xu = x;
        if  x(ku) == 0,  xu(ku) = opts(5)^2;
        else,            xu(ku) = x(ku) + opts(5)*abs(x(ku)); end
        [stop Fu fu] = checkfJ(fun,xu,varargin{:});  info(7) = info(7)+1;
        if  ~stop
          hu = xu - x;
          B = B + ((fu - f - B*hu)/norm(hu)^2) * hu';
        end
      end
      B = B + ((fn - f - B*h)/norm(h)^2) * h'; 
      k = k + 1;
      if  Trace
        X(:,k) = xnew;   perf(:,k) = [Fn norm(B'*fn,inf) mu]'; end
      dL = (h'*(mu*h - g))/2;  dF = F - Fn;
      if  (dL > 0) & (dF > 0)               % Update x and modify mu
        x = xnew;   F = Fn;  f = fn;
        mu = mu * max(1/3, 1 - (2*dF/dL - 1)^3);   nu = 2;
      else  % Same  x, increase  mu
        mu = mu*nu;  nu = 2*nu; 
      end 
      if  k > kmax,  stop = 3; 
      else
        g = B'*f;  ng = norm(g,inf);      A = B'*B;
        if  isinf(ng) | isinf(norm(A(:),inf)),  stop = -5; end
      end
    end  
  end
end
%  Set return values
if  Trace
  X = X(:,1:k);   perf = perf(:,1:k);
else,  X = x;  end
if  stop < 0,  tau = NaN;  else,  tau = mu/max(diag(A)); end
info(1:6) = [F  ng  nh  tau  k-1  stop];