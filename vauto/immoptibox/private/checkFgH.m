function  [err, F,g,H] = checkFgH(fun,x,varargin)
%CHECKFG  Check Matlab function which is called by a 
% general optimization function

% Version 04.04.03.  hbn@imm.dtu.dk

err = 0; 
[F g H] = feval(fun,x,varargin{:});
sF = size(F);   sg = size(g);   sH = size(H);
if  any(sF ~= 1) | ~isreal(F) | any(isnan(F(:))) | any(isinf(F(:)))
  err = -2;  return, end
if  ~isreal(g) | any(isnan(g(:))) | any(isinf(g(:)))
  err = -3;  return, end
if  min(sg) ~= 1 | max(sg) ~= length(x)
  err = -4;  return, end
if  ~isreal(H) | any(isnan(H(:))) | any(isinf(H(:)))
  err = -3;  return, end
if  any(sH ~= length(x))
  err = -4;  return, end
E = H - H';
if  norm(E(:)) > 10*eps*norm(H(:))
  err = -5; 
end
H = (H + H')/2;