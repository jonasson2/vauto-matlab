function  [X, info, perf] = dampnewton(fun, x0, opts, varargin)
%DAMPNEWTON  Damped Newton method for unconstrained optimization:
% Find  xm = argmin{F(x)} , where  x  is an n-vector and the scalar
% function  F  with gradient  g  (with elements  g(i) = DF/Dx_i )
% and Hessian  H  (with elements  H(j,k) = D^2F/(Dx_j Dx_k) )
% must be given by a MATLAB function with with declaration
%            function  [F, g, H] = fun(x,p1,p2,...)
% p1,p2,... are parameters of the function.
%  
% Call
%    [X, info] = dampnewton(fun, x0)
%    [X, info] = dampnewton(fun, x0, opts, p1,p2,...)
%    [X, info, perf] = dampnewton(.....)
%
% Input parameters
% fun  :  Handle to the function.
% x0   :  Starting guess for  xm .
% opts :  Vector with four elements.  
%         opts(1)  used in starting value for damping parameter: 
%             mu = opts(1) * max(abs(H0(i,i))), where H0 is the 
%             Hessian of F at x0.
%         opts(2:4)  used in stopping criteria:
%             ||F'||inf <= opts(2)                     or 
%             ||dx||2 <= opts(3)*(opts(3) + ||x||2)    or
%             no. of iteration steps exceeds  opts(4) .
%         Default  opts = [1e-6 1e-4 1e-8 100]
%         If the input opts has less than 4 elements, it is
%         augmented by the default values.
% p1,p2,..  are passed directly to the function FUN .    
%
% Output parameters
% X    :  If  perf  is present, then array, holding the iterates
%         columnwise.  Otherwise, computed solution vector.
% info :  Performance information, vector with 6 elements:
%         info(1:4) = final values of 
%             [F(x)  ||F'||inf  ||dx||2  mu/max(H(i,i))] .
%         info(5) = no. of iteration steps
%         info(6) = 1 : Stopped by small gradient
%                   2 :  Stopped by small x-step
%                   3 :  No. of iteration steps exceeded 
%                  -1 :  x is not a real valued vector
%                  -2 :  F is not a real valued scalar
%                  -3 :  g is not a real valued vector or
%                        H is not a real valued matrix
%                  -4 :  Dimension mismatch in x, g, H
%                  -5 :  H is not symmetric
% perf :  Array, holding 
%         perf(1,:) = values of  F(x)
%         perf(2,:) = values of  || F'(x) ||inf
%         perf(3,:) = mu-values.

% Version 04.09.07.  hbn(a)imm.dtu.dk

% Check parameters and function call
if  nargin < 3,  opts = []; end
opts = checkopts(opts, [1e-6 1e-4 1e-8 100]); 

if  nargin < 2,  stop = -1;
else
  [stop x n] = checkx(x0);   
  if  ~stop,  [stop F g H] = checkFgH(fun,x0,varargin{:}); end
end
if  stop
  X = x0;  perf = [];  info = [repmat(NaN,1,4) 0  stop];
  return
end

%  Finish initialization
% --------------  Modified 04.09.07
mu = opts(1) * norm(H(:),inf);  
if  mu == 0  % zero initial Hessian.  Steepest descent direction
  mu = norm(g)/max(norm(x), sqrt(eps));
end
% --------------

ng = norm(g,inf);   kmax = opts(4);
Trace = nargout > 2;
if  Trace
  X = repmat(x,1,kmax+1);
  perf = repmat([F; ng; mu],1,kmax+1);
end 
k = 1;   nu = 2;   nh = 0;   stop = 0;

% Iterate
while   ~stop 
  if  ng <= opts(2),  stop = 1;  
  else
    [h mu] = geth(H,g,mu);
    nh = norm(h);   nx = opts(3) + norm(x);
    if  nh <= opts(3)*nx,  stop = 2; end 
  end 
  if  ~stop
    xnew = x + h;   h = xnew - x;   
    [stop Fn gn Hn] = checkFgH(fun,xnew,varargin{:}); 
    if  ~stop      
      k = k + 1;  dF = F - Fn;   accept = 0;
      if  Trace
        X(:,k) = xnew;   perf(:,k) = [Fn norm(gn,inf) mu]'; end
      if  dF > 0,  accept = 1;  dL = (h'*(mu*h - g))/2;
      elseif  Fn <= F + abs(F)*(1 + 100*eps)  % Try gradient
        dF = (g + gn)'*(g - gn);
        if  dF > 0,  accept = 2; end
      end        
      if  accept                     % Update x and modify mu
        if  accept == 1 & dL > 0
          mu = mu * max(1/3, 1 - (2*dF/dL - 1)^3);  nu = 2;
        else
          mu = mu*nu;  nu = 2*nu;  
        end
        x = xnew;   F = Fn;  g = gn;  H = Hn;   ng = norm(g,inf);
      else                                  % Same  x, increase  mu
        mu = mu*nu;  nu = 2*nu; 
      end
      if  k > kmax, stop = 3;  end
    end
  end   
end
%  Set return values
if  Trace
  X = X(:,1:k);   perf = perf(:,1:k);
else,  X = x;  end
if  stop < 0,  F = NaN;  ng = NaN; end
info = [F  ng  nh  mu/max(abs(diag(H)))  k-1  stop];