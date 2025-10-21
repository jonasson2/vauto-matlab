function  [X, info, perf] = dogleg(fun, x0, opts, varargin)
%DOGLEG Powell's dog-leg method for least squares.
% Find  xm = argmin{F(x)} , where  x  is an n-vector and
% F(x) = .5 * sum(f_i(x)^2) .
% The functions  f_i(x) (i=1,...,m) and the Jacobian matrix  J(x)  
% (with elements  J(i,j) = Df_i/Dx_j ) must be given by a MATLAB
% function with declaration
%            function  [f, J] = fun(x, p1,p2,...)
% p1,p2,... are parameters of the function.  In connection with 
% nonlinear data fitting they may be arrays with coordinates of 
% the data points.
%  
% Call
%    [X, info] = dogleg(fun, x0)
%    [X, info] = dogleg(fun, x0, opts, p1,p2,...)
%    [X, info, perf] = dogleg(.....)
%
% Input parameters
% fun  :  Handle to the function.
% x0   :  Starting guess for  xm .
% opts :  Vector with five elements:
%         opts(1)    initial trust region radius Delta
%         opts(2:5)  used in stopping criteria:
%             ||F'||inf <= opts(2)                     or 
%             ||dx||2 <= opts(3)*(opts(3) + ||x||2)    or
%             ||f||inf <= opts(4)                      or
%             no. of iteration steps exceeds  opts(5) . 
%         Default  opts = [0.1(1+||x0||) 1e-4 1e-8 1e-6 100]
%         If the input opts has less than 5 elements, it is
%         augmented by the default values.
% p1,p2,..  are passed dirctly to the function FUN .    
%
% Output parameters
% X    :  If  perf  is present, then array, holding the iterates
%         columnwise.  Otherwise, computed solution vector.
% info :  Performance information, vector with 6 elements:
%         info(1:4) = final values of 
%             [F(x)  ||F'||inf  ||dx||2  Delta] 
%         info(5) = no. of iteration steps
%         info(6) = 1 :  Stopped by small gradient
%                   2 :  Stopped by small x-step
%                   3 :  Stopped by small f-vector
%                   4 :  No. of iteration steps exceeded
%                  -1 :  x is not a real valued vector
%                  -2 :  f is not a real valued column vector
%                  -3 :  J is not a real valued matrix
%                  -4 :  Dimension mismatch in x, f, J
%                  -5 :  Overflow during computation 
% perf :  Array, holding 
%         perf(1,:) = values of  F(x)
%         perf(2,:) = values of  || F'(x) ||inf
%         perf(3,:) = Delta-values.

% Version 04.09.08.  hbn(a)imm.dtu.dk

% Check parameters and function call
if  nargin < 2,  stop = -1;
else
  [stop x n] = checkx(x0);   
  if  ~stop,  [stop F f J] = checkfJ(fun,x0,varargin{:}); end
end
if  ~stop
  g = -(J'*f);   ngi = norm(g,inf);
  if  isinf(ngi),  stop = -5; end
else
  F = NaN;  ngi = NaN;
end
if  stop
  X = x0;  perf = [];  info = [F  ngi  0  opts(1)  0  stop];
  return
end
  
% Finish initialization
if  nargin < 3,  opts = []; end
opts = checkopts(opts, [.1*(1 + norm(x)) 1e-4 1e-8 1e-6 100]); 
Delta = opts(1);    kmax = opts(5);
Trace = nargout > 2;
if  Trace
  X = repmat(x,1,kmax+1);
  perf = repmat([F; ngi; Delta],1,kmax+1);
end 
k = 1;   nu = 2;   nx = norm(x);  fact = 1;  stop = 0;
reduce = 0;   nstep = 0;
% Iterate
while   ~stop 
  if      isinf(ngi),                       stop = -5;
  elseif  ngi <= opts(2),                   stop = 1;  
  elseif  Delta <= opts(3)*(opts(3) + nx),  stop = 2;  
  elseif  norm(f,inf) <= opts(4),           stop = 3; 
  elseif  k >= kmax,                        stop = 4;
  else
    if  fact  % Factorize and ompute hGN
      [Q R] = qr(J,0);  [U S V] = svd(R);
      s = diag(S);  i = find(s > 100*eps*s(1));
      Qf = -(f'*Q)';   UQf = U(:,i)'*Qf;
      hGN = V(:,i)*(UQf ./ s(i));  
      nhGN = norm(hGN);  fact = 0;
    end
    if  nhGN > Delta  % Include gradient
      nh = Delta;
      ng = norm(g);  alpha = (ng / norm(J*g))^2;
      gn = alpha*g;  ngn = alpha*ng;
      if  ngn >= Delta
        h = (Delta/ngn) * gn;
        dLpre = Delta*(2*ngn - Delta)/(2*alpha);
      else  %  Dog leg
        b = hGN - gn;  bb = b'*b;  gb = gn'*b;
        c = (Delta + ngn)*(Delta - ngn);
        if  gb > 0
          beta = c / (gb + sqrt(gb^2 + c * bb));
        else
          beta = (sqrt(gb^2 + c * bb) - gb)/bb;
        end
        h = gn + beta*b;
        dLpre = .5*alpha*(1 - beta)^2*ng^2 + beta*(2-beta)*F;
      end
    else
      h = hGN;  nh = nhGN;
      dLpre = F;  
      if  nh <= opts(3)*(opts(3) + nx),  stop = 2; end
    end
  end
  if  ~stop
    xnew = x + h;   h = xnew - x;   nstep = norm(h);
    dL = F - .5*norm(f + J*h)^2;
    [stop Fn fn Jn] = checkfJ(fun,xnew,varargin{:});
    if  ~stop
      dF = F - Fn;
      if  dF > 0 & dL > 0                        % Update x 
        x = xnew;   F = Fn;  J = Jn;  f = fn;   fact = 1;
        g = -(J'*f);   ngi = norm(g,inf);
        rho = dL / dF;
      else,  rho = -1; end
      % Update Delta
      if  abs(rho-1) < .2 & nh > Delta/3 & reduce <= 0
        Delta = 3*Delta;   nu = 2;  reduce = 0;
      elseif  rho < .25
        Delta = Delta/nu;
        nu = 2*nu;  reduce = 2;
      else
        reduce = reduce - 1;
      end  
      k = k + 1;
      if  Trace,  X(:,k) = xnew;   perf(:,k) = [Fn ngi Delta]'; end
      if  k > kmax,  stop = 3; end 
    end
  end   
end
%  Set return values
if  Trace
  X = X(:,1:k);   perf = perf(:,1:k);
else,  X = x;  end
if  stop < 0,  F = NaN;  ngi = NaN; end
info = [F  ngi  nstep  Delta  k-1  stop];