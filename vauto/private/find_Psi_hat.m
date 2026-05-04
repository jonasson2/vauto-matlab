function Psi_hat = find_Psi_hat(Psi, Sig)
  LSig = chol(Sig, 'lower');
  r = size(Sig, 1);
  h = size(Psi, 1)/r;
  I = 1:r;
  Psi_hat = zeros(r*h, r*h);
  for k = 1:h
    Psi_hat(:,I) = Psi(:,I)*LSig;
    I = I+r;
  end
end
