function [Q, R] = weighted_qr(K, W)
%WEIGHTED_QR  Thin QR factorization, orthonormal in the W-inner product.
%   [Q,R] = WEIGHTED_QR(K, W) returns K = Q*R with Q'*W*Q = I_r, where K is
%   m x r (m >= r), Q is m x r and R is r x r upper triangular. W is an SPD
%   weight (vector of diagonal entries, diagonal matrix, or full SPD
%   matrix, see WEIGHT_FACTOR).
%
%   This is QR_M (W = M, mass matrix) and QR_Xi (W = Xi, quadrature
%   weights) of the paper: scale, plain QR, unscale.

[apply, unapply] = weight_factor(W);
[Q, R] = qr(apply(K), 0);
Q = unapply(Q);
end