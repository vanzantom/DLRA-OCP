function [U, S, V] = weighted_svd(Y, r, M, Xi)
%WEIGHTED_SVD  Best rank-r approximation in the (M,Xi)-weighted norm.
%   [U,S,V] = WEIGHTED_SVD(Y, r, M, Xi) returns the truncated SVD of the
%   m x n matrix Y with respect to the inner product <X,Y> = tr(X'*M*Y*Xi):
%       Y ~ U*S*V',  U'*M*U = I_r,  V'*Xi*V = I_r,  S = diag(s_1..s_r).
%   M and Xi are SPD weights (vectors, diagonal or full matrices, see
%   WEIGHT_FACTOR). If Y has fewer than r singular values, r is reduced.
%
%   Algorithm: Yh = C_M * Y * C_Xi' (C'*C = weight), plain SVD, unscale.

[applyM, unapplyM]   = weight_factor(M);
[applyXi, unapplyXi] = weight_factor(Xi);

Yh = applyXi(applyM(Y)')';           % = C_M * Y * C_Xi'
[Uh, Sh, Vh] = svd(full(Yh), 'econ');
r = min(r, size(Sh, 1));
U = unapplyM(Uh(:, 1:r));
V = unapplyXi(Vh(:, 1:r));
S = Sh(1:r, 1:r);
end