function  xx = delknot(x,k)
% Reduced knot set when knot no k is deleted

% Version 04.04.28.  hbn@imm.dtu.dk

h = diff(x);   n = length(h);
xx = x([1:k-1 k+1:end]);   hh = h(k-1) + h(k);
% Adjust adjoining intervals
if  k > 2,  xx(k-1) = xx(k-1) + .1*min(h(k-2),hh); end
if  k < n-1,  xx(k) = xx(k) - .1*min(h(k+1),hh); end