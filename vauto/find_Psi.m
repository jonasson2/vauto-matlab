function Psi = find_Psi(A, B)
  % Get dimensions and prepare
  [r, p, q] = getdims(A, B);
  h = max(p, q);
  if h == 0, Psi = []; return, end
  Aflp = flipmat(A);
  Psi = zeros(h*r, r);

  % Compute first block column
  J = 1:r;
  k = p*r + 1;
  l = 0;
  Psi(J,:) = eye(r);
  for j = 2:h
    J = J + r;
    k = k - r;
    l = l + r;
    Psi(J,:) = Aflp(:, k:end)*Psi(1:l,:) + B(J - r);
  end

  % Copy to remaining lower Psi blocks
  J = 1:r;
  k = h*r;
  l = 1;
  for j = 2:h
    J = J + r;
    k = k - r;
    l = l + r;
    Psi(l:end, J) = Psi(1:k, 1:r);
  end
  
  % Copy to upper triangle (use that diagonal blocks are I).
  Psi = Psi + tril(Psi,-1)';
end