function G = find_G(A, Bflp, Psi, Sig)
  % Get dimensions
  [r,p,q] = getdims(A, Bflp);
  h = max(p,q);

  % Copy last p blocks from bottom block row of Psi
  I = h*r-r+1 : h*r;
  k = h*r - p*r + 1;
  Psi_tilde = Psi(I, k:end);
  
  % Post-multiply with Sig
  K = 1:r;
  for k = 1:p
    Psi_tilde(:,K) = Psi_tilde(:,K)*Sig;
    K = K + r;
  end
  
  % Determine G = [G0...G{p-1}]
  J = 1:r;
  G = zeros(r, r*p);
  G(:,J) = Sig + Bflp*Psi_tilde';
  k = 1;
  for j = 2:min(p,q)
    J = J + r;
    G(:,J) = Bflp(:,k:end)*Psi_tilde(:,k:end)';
    k = k + r;
  end
end