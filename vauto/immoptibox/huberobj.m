function  [F,S,f,J,g] = huberobj(fun,x,gamma,Htype,varargin)
%HUBEROBJ  Value and gradient of the Huber objective function
% for the vector function given by
%            function  [f, J] = fun(x,p1,p2,...)
% p1,p2,... are parameters of the function and J is the Jacobean.
% If nargout < 4, then the function only needs to return f, and the
% gradient is not computed.
%  
% Call
%    [F,S,f] = huberobj(fun,x,gamma)
%    [F,S,f] = huberobj(fun,x,gamma,Htype, p1,p2,...)
%    [F,S,f,J,g] = huberobj(.....)
%
% Input parameters
% fun   :  Handle to the function.
% x     :  n-vector, argument of fun
% gamma :  Huber threshold.
% Htype :  Choice of Huber function,
%          1 : one-sided,  f_i > 0  are neglected,
%          2 : one-sided,  all f_i <= gamma are active,
%          otherwise,  all abs(f_i) <= gamma  are active (default).
% p1,p2,..  are passed directly to the function FUN .    
%
% Output parameters
% F   :  Huber objective function.
% S   :  Struct with the Huber active set.  Fields
%        S.s  Huber sign vector,
%        S.A  active set (indices of small components in f),
%        S.N  indices of large components,
%        S.L  = [length(S.A), length(S.N)].
% f,J :  Output from the evaluation of fun
% g   :  Gradient,  g = F'(x).   

% Version 04.09.05.  hbn(a)imm.dtu.dk

% Evaluate the function and set rounding error threshold
if  nargout < 4
  f = feval(fun,x, varargin{:});
  tau = 10*eps*norm(f,inf);
else
  [f J] = feval(fun,x, varargin{:});
  tau = 10*eps*max(norm(f,inf), norm(J,inf)*norm(x,inf));
end
if  gamma < tau
  error(sprintf('gamma must be at least %9.2e',tau))
end
thr = gamma + tau;

% Huber sign vector
s = sign(f);
if  nargin < 4 | abs(Htype - 1.5) ~= 0.5,  Htype = 3; end
if  Htype == 1  % Neglect positive contributions
  w = find(-thr <= f & f <= tau);    s(w) = 0;
  p = find(s > 0);   s(p) = 0;
elseif  Htype == 2  % Negative contributions are active
  w = find(f <= thr);   s(w) = 0;
else  % simple Huber function
  w = find(abs(f) <= thr);    s(w) = 0; 
end

% Active set
N = find(s);
S = struct('s',s', 'A',w', 'N',N', 'L',[length(w)  length(N)]);

% Compute function
fA = f(S.A);  fN = f(S.N);  
F = norm(fA)^2/(2*gamma) + norm(fN,1) - .5*gamma*S.L(2);

if  nargout > 3  % gradient with check of rounding errors
  g = J(S.A,:)'*fA/gamma + J(S.N,:)'*s(S.N);
  thr = 10*eps*(norm(J(S.A,:),1)*norm(fA,inf)/gamma + norm(J(S.N,:),1));
  if  norm(g,inf) <= thr,  g = zeros(size(g)); end  
end  % gradient