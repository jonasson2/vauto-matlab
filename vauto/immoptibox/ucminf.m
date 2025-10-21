function  [X, info, perf, D] = ucminf(fun, x0, opts, D0, varargin)
%UCMINF  BFGS method for unconstrained nonlinear optimization:
% Find  xm = argmin{F(x)} , where  x  is an n-vector and the scalar
% function  F  with gradient  g  (with elements  g(i) = DF/Dx_i )
% must be given by a MATLAB function with with declaration
%            function  [F, g] = fun(x,p1,p2,...)
% p1,p2,... are parameters of the function.
% 
% Call
%    [X, info] = ucminf(fun, x0)
%    [X, info] = ucminf(fun, x0, opts)
%    [X, info] = ucminf(fun, x0, opts, D0, p1,p2,...)
%    [X, info, perf] = ucminf(.....)
%    [X, info, perf, D] = ucminf(.....)
%
% Input parameters
% fun  :  Handle to the function.
% x0   :  Starting guess for  xm .
% opts :  Vector with four elements.  
%         opts(1) :  Expected length of initial step
%         opts(2:4)  used in stopping criteria:
%             ||g||_inf <= opts(2)                     or 
%             ||dx||_2 <= opts(3)*(opts(3) + ||x||_2)  or
%             no. of function evaluations exceeds  opts(4) . 
%         Default  opts = [1  1e-4  1e-8  100]
%         If the input opts has less than 4 elements, it is
%         augmented by the default values.
% D0   :  If present, then approximate inverse Hessian at  x0 .
%         Otherwise, D0 := I.
% p1,p2,..  are passed dirctly to the function FUN .
%
% Output parameters
% X    :  If  perf  is present, then array, holding the iterates
%         columnwise.  Otherwise, computed solution vector.
% info :  Performance information, vector with 6 elements:
%         info(1:3) = final values of [F(x)  ||g||_inf  ||dx||_2] 
%         info(4:5) = no. of iteration steps and evaluations of (F,g)
%         info(6) = 0 :  Line search failed
%                   1 :  Stopped by small gradient
%                   2 :  Stopped by small x-step
%                   3 :  No. of iteration steps exceeded 
%                   4 :  Stopped by zero step.
%                  -1 :  x is not a real valued vector
%                  -2 :  F is not a real valued scalar
%                  -3 :  g is not a real valued vector 
%                  -4 :  x and g have different lengths
%                  -6 :  D0 is not a symmetric, pos. def. n*n-matrix
% perf :  Array, holding 
%         perf(1:2,:) = values of  F(x) and ||g||_inf
%         perf(3,:)   = trust region radius.
%         perf(4:6,:) = Line search info, cf LINESEARCH:  values of  
%                       am, U'(am), no. fct. evals.
% D   :   Array holding the approximate inverse Hessian at 
%         the computed minimizer.

% Version 06.11.08.  hbn(a)imm.dtu.dk

% Initial check
if  nargin < 2,  error('Too few input parameters'), end

% Check OPTS
if  nargin < 3,  opts = []; end
opts = checkopts(opts, [1 1e-4 1e-8 100]); 

% Check parameters and function call
[stop x n] = checkx(x0);   
if  ~stop
  [stop F g] = checkfg(fun,x0,varargin{:}); 
  if  ~stop  % Initial inverse Hessian
    if  nargin > 3 & ~isempty(D0)
      [stop D] = checkD(n,D0);  fst = 0;
    else,   D = eye(n);  fst = 1;    end
  end
else,  F = NaN;  end

if  stop
  X = x0;  perf = [];  D = [];  info = [F(1)  NaN  0  0  1  stop];
  return
end

%  Finish initialization
k = 1;   kmax = opts(4);   neval = 1;   ng = norm(g,inf);
Delta = opts(1);
Trace = nargout > 2;
if  Trace
  X = repmat(x(:),1,kmax+1);
  perf = repmat([F; ng; Delta; zeros(3,1)],1,kmax+1);
end
if  ng <= opts(2),  stop = 1;  nh = 0;
else
  h = zeros(size(x));  nh = 0;
  ngs = repmat(ng,1,3);
  lsopts = [1 .05 .99 5 2];  % LINESEARCH options
end

more = 1;
while  ~stop & more
  %  Previous values
  xp = x;   gp = g;   Fp = F;   nx = norm(x);
  ngs = [ngs(2:3) ng];
  h = D*(-g(:));   nh = norm(h);   red = 0; 
  if  nh <= opts(3)*(opts(3) + nx),  stop = 2;  
  else
    if  fst | nh > Delta  % Scale to ||h|| = Delta
      h = (Delta / nh) * h;   nh = Delta;   
      fst = 0;  red = 1;
    end
    %  Line search
    [x F g  linfo] = linesearch(fun,x,F,g, h, lsopts, varargin{:}); 
    neval = neval + linfo(3);  
    if  linfo(1) <= 0
      stop = linfo(1);  more = 0;  % something wrong
    else
      k = k+1;
      if  linfo(1) < 1  % Reduce Delta
        Delta = .35 * Delta;
      elseif   red & (linfo(2) > .7)  % Increase Delta
        Delta = 3*Delta;      
      end 
      %  Update ||g||
      ng = norm(g,inf);
      if  Trace
        X(:,k) = x(:); 
        perf(:,k) = [F; ng; Delta; linfo(1); dot(h,g); linfo(3)]; 
      end
      h = x - xp;   nh = norm(h);
      if  nh == 0,
        stop = 4; 
      else
        y = g - gp;   yh = dot(y,h);
        if  yh > sqrt(eps) * nh * norm(y)
          %  Update  D
          v = D*y(:);   yv = dot(y,v);
          a = (1 + yv/yh)/yh;   w = (a/2)*h(:) - v/yh;
          D = D + w*h' + h*w';
        end  % update D
        %  Check stopping criteria
        thrx = opts(3)*(opts(3) + norm(x));
        if      ng <= opts(2),              stop = 1;
        elseif  nh <= thrx,                 stop = 2;
        elseif  neval >= kmax,              stop = 3; 
        else,   Delta = max(Delta, 2*thrx);  end
      end  
    end  % Nonzero h
  end % nofail
end  % iteration

%  Set return values
if  Trace
  X = X(:,1:k);   perf = perf(:,1:k);
else,  X = x;  end
info = [F  ng  nh  k-1  neval  stop];

% ==========  auxiliary function  =================================

function  [err, D] = checkD(n,D0)
% Check given inverse Hessian
D = D0;   sD = size(D);  err = 0;
if  any(sD - n),  err = -6;  return, end
% Check symmetry
dD = D - D';   ndD = norm(dD(:),inf);
if  ndD > 10*eps*norm(D(:),inf),  err = -6;  return, end
if  ndD,  D = (D + D')/2; end  % Symmetrize      
[R p] = chol(D);
if  p,  err = -6;  end