function  [X, info, perf] = nlshybrid(fun, x0, opts, varargin)
%NLSHYBRID  Levenberg-Marquardt and Qausi-Newton hybrid method for 
% nonlinear least squares.
% Find  xm = argmin{F(x)} , where  x = [x_1, ..., x_n]  and
% F(x) = .5 * sum(f_i(x)^2) .
% The functions  f_i(x) (i=1,...,m) and the Jacobian matrix  J(x)  
% (with elements  J(i,j) = Df_i/Dx_j ) must be given by a MATLAB
% function with declaration
%            function  [f, J] = fun(x,p1,p2,...)
% p1,p2,... are parameters of the function.  In connection with 
% nonlinear data fitting they may be arrays with coordinates of 
% the data points.
% 
% Call
%    [X, info] = nlshybrid(fun, x0)
%    [X, info] = nlshybrid(fun, x0, opts, p1,p2,...)
%    [X, info, perf] = nlshybrid(......)
%
% Input parameters
% fun  :  Handle to the function.
% x0   :  Starting guess for  xm .
% opts :  Vector with four elements.  
%         opts(1)  used in starting value for Marquardt parameter: 
%             mu = opts(1) * max(A0(i,i))  with  A0 = J(x0)'*J(x0)
%         opts(2:4)  used in stopping criteria:
%             ||F'||inf <= opts(2)                     or 
%             ||dx||2 <= opts(3)*(opts(3) + ||x||2)    or
%             no. of iteration steps exceeds  opts(4) .
%         Default  opts = [1e-3 1e-4 1e-8 100]
%         If the input opts has less than 4 elements, it is
%         augmented by the default values.
% p1,p2,..  are passed directly to the function FUN .   
%
% Output parameters
% X    :  If  perf  is present, then array, holding the iterates
%         columnwise.  Otherwise, computed solution vector.
% info :  Performance information, vector with 6 elements:
%         info(1:4) = final values of 
%             [F(x)  ||F'||inf  ||dx||2  mu/max(A(i,i))] ,
%           where  A = J(x)'*J(x) .
%         info(5) = no. of iteration steps
%         info(6) = 1 : Stopped by small gradient
%                   2 :  Stopped by small x-step
%                   3 :  No. of iteration steps exceeded 
%                  -1 :  x is not a real valued vector
%                  -2 :  f is not a real valued column vector
%                  -3 :  J is not a real valued matrix
%                  -4 :  Dimension mismatch in x, f, J
%                  -5 :  Overflow during computation 
% perf :  Array, holding 
%         perf(1,:) = values of  F(x)
%         perf(2,:) = values of  || F'(x) ||inf
%         perf(3,:) = mu-values
%         perf(4,:) = method (1 for L-M, 2 for Q-N)

% Version 04.03.19.  hbn(a)imm.dtu.dk

% Check parameters and function call
if  nargin < 3,  opts = []; end
opts = checkopts(opts, [1e-3 1e-4 1e-8 100]); 
if  nargin < 2,  stop = -1;
else
  [stop x n] = checkx(x0);   
  if  ~stop,  [stop F f J] = checkfJ(fun,x0,varargin{:}); end
end
if  ~stop
  g = J'*f;   ng = norm(g,inf);  A = J'*J;
  if  isinf(ng) | isinf(norm(A(:),inf)),  stop = -5; end
else
  F = NaN;  ng = NaN;
end
if  stop
  X = x0;  perf = [];  info = [F  ng  0  opts(1)  0  stop];
  return
end

%  Finish initialization
mu = opts(1) * max(diag(A));    kmax = opts(4);
Trace = nargout > 2;
if  Trace
  X = repmat(x,1,kmax+1);
  perf = repmat([F; ng; mu; 1],1,kmax+1);
end 

% Struct with information about a point
pt = struct('x',x, 'f',f, 'J',J, 'g',g, 'F',F, ...
  'ng',[ng  repmat(NaN,1,3)], ...  % four recent values of ||g||
  'muD',[mu 2 0], ...              % mu, nu, Delta
  'aux',[size(J) 0 1]);            % [m  n  count  method]

% Initial Hessian approximation
B = eye(pt.aux(2));   

% Iterate
k = 1;   stop = 0;   kmax = opts(4);
while  ~stop
  if  pt.ng(1) <= opts(2),  stop = 1;
  else   
    if  pt.aux(4) == 1
      [ptn  stop  better nh] = LMstep(fun,pt,opts,varargin{:});
    else
      [ptn  stop  better nh] = QNstep(fun,pt,opts,B,varargin{:});
    end
    if  stop,  pt = ptn;   
    else     % Continue
      if  ptn.F < 1.5*pt.F   %  Update  B
        h = ptn.x - pt.x;
        if  norm(h)
          y = ptn.J'*(ptn.J*h) + (ptn.J - pt.J)'*ptn.f;  
          hy = dot(h,y);
          if  hy > 0
            v = B*h;
            B = B + y*(y/hy)' - v*(v/dot(h,v))';
          end
        end
      end  % update B
      if  better,  pt = ptn;   % update x
      else,        pt.muD = ptn.muD;  end
    end 
    k = k + 1;
    if  Trace 
      X(:,k) = pt.x;
      perf(:,k) = [pt.F pt.ng(1) pt.muD(1) pt.aux(4)]';
    end
  end
  if  k > kmax,  stop = 3; end 
end
% Set return values
if  Trace
  X = X(:,1:k);   perf = perf(:,1:k);
else,  X = pt.x;  end 
newopts1 = pt.muD(1)/max(sum(pt.J .* pt.J));
info = [pt.F  pt.ng(1)  nh  newopts1  k-1  stop];

% ==========  auxiliary functions  =================================

function  [ptn, stop, better, nh] = LMstep(fun,pt,opts,varargin)
% Levenberg-Marquardt step
ptn = pt;   stop = 0;  better = 0;  method = 1;  
mu = pt.muD(1);   count = pt.aux(3); 
% disp([mu count])
h = geth(pt.J'*pt.J, pt.g, mu);
nh = norm(h);   nx = opts(3) + norm(pt.x); 
if      nh <= opts(3)*nx,  stop = 2;
else
  ptn.x = pt.x + h;   
  [stop ptn.F ptn.f ptn.J] = checkfJ(fun, ptn.x, varargin{:}); 
  if  ~stop
    ptn.g =  ptn.J'*ptn.f;   
    ptn.ng = [norm(ptn.g,inf) pt.ng(1:end-1)];
    if  ptn.ng(1) <= opts(2),  stop = 1;  % small gradient
    else
      h = ptn.x - pt.x;   dL = (h'*(mu*h - pt.g))/2;  % predicted decrease
      dF = pt.F - ptn.F;                  % actual decrease  
      better = dF > 0;
      if  (dL > 0) & (dF > 0)
        if  ptn.ng(1) < .02*ptn.F
          count = count + 1;  % relatively small gradient
        else
          count = 0;  
        end
        if  count >= 3,  method = 2;  end  % switch method
        ptn.muD(1:2) = [mu * max(1/3, (1 - (2*dF/dL - 1)^3))  2];
      else
        ptn.muD(1:2) = pt.muD(1:2).*[pt.muD(2) 2];
        count = 0;
      end
      q = ptn.ng ./ pt.ng;
      if  ~any(isnan(q)) & all(q < 1 & q > .8) & norm(diff(q),inf) < 1e-3
        method = 2;   % slow, linear convergence
      end
    end
  end   
end
if  method == 2  % set trust region radius
  ptn.muD(3) = max(1.5*opts(3)*nx, nh/5);
%   D_thr = [ptn.muD(3)  opts(3)*nx]
%   LMx = ptn.x
%   opts
end
ptn.aux(3:4) = [count method];   
  
function  [ptn, stop, better, nh] = QNstep(fun,pt,opts,B,varargin)
% Quasi-Newton step  
ptn = pt;   stop = 0;   better = 0;   method = 2;
nx = opts(3) + norm(pt.x); 
if  pt.muD(3) <= opts(3)*nx
  stop = 2;  
else
  h = geth(B,pt.g, eps*norm(diag(B),inf));
  nh = norm(h);  
  if   nh <= opts(3)*nx,  stop = 2; end
end
if  stop,  return,  end

% If necessary, shrink to trust region
if  nh > pt.muD(3)
  h = (pt.muD(3)/nh)*h;   nh = pt.muD(3);
end
ptn.x = pt.x + h;
[stop ptn.F ptn.f ptn.J] = checkfJ(fun, ptn.x, varargin{:}); 
if  ~stop
  ptn.g = ptn.J'*ptn.f;   
  ptn.ng = [norm(ptn.g,inf) pt.ng(1:end-1)];
  better = (ptn.F < pt.F) | ...
    (ptn.F <= (1+sqrt(eps))*pt.F & ptn.ng(1) < pt.ng(1));
  if      ptn.ng(1) <= opts(2),  stop = 1;     % small gradient
  elseif  ptn.ng(1) >= pt.ng(1),  method = 1;  % no decrease in gradient
  else  % Update trust region
    h = ptn.x - pt.x;   nh = norm(h);
    dL = -dot(h, pt.g)/2;     % predicted decrease
    dF = pt.F - ptn.F;        % actual decrease  
    if  dL,  rho = dF/dL; else,  rho = .5; end
    if      rho < .25,  ptn.muD(3) = ptn.muD(3)/2;
    elseif  rho > .75,  ptn.muD(3) = max(ptn.muD(3), 3*nh); end
  end
end 
ptn.aux(3:4) = [0  method];