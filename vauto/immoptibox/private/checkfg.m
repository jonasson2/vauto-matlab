function  [err, f,g] = checkfg(fun,x,varargin)
%CHECKFG  Check Matlab function which is called by a 
% general optimization function

% Version 04.01.26.  hbn@imm.dtu.dk

err = 0; 
[f g] = feval(fun,x,varargin{:});
sf = size(f);   sg = size(g);
if  any(sf ~= 1) | ~isreal(f) | any(isnan(f(:))) | any(isinf(f(:)))
  err = -2;  return, end
if  ~isreal(g) | any(isnan(g(:))) | any(isinf(g(:)))
  err = -3;  return, end
if  min(sg) ~= 1 | max(sg) ~= length(x)
  err = -4;  return, end