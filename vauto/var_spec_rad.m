%   VAR_SPEC_RAD(A) returns the spectral radius of the companion matrix, Ac, of
%   the VAR process
%                  x(t) = A1·x(t-1) + ... + Ap·x(t-p) + eps(t)
%
%   where each Ai is r×r and A = [A1 A2...Ap]. The companion matrix is 
%
%       Ac = A1  A2  ... Ap-1  Ap
%             I
%                 I
%                    ...
%                          I    O
%
%   and its specral radius, rho = max(eig(Ac)) measures the asymptotic growth
%   rate of the process excluding shocks, the process being stationary when rho
%   < 1. 
% 
% NOTE: The stability of a VARMA process with coefficents A1, A2... and B1,
% B2... depends only on the Ai matrices.

function rho = var_spec_rad(A)
  r = size(A, 1);
  p = length(A)/r;
  k = r*(p-1);
  Ac = [A; eye(k) zeros(k,r)];
  rho = max(abs(eig(Ac)));
end