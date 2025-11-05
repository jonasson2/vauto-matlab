function Aflp = flipmat(A)
  % Flips [A1...An] to [An...A1]
  r = size(A, 1);
  n = size(A, 2)/r;
  Aflp reshape(flipdim(reshape(A,r,r,n), 3), r, r*n);
end