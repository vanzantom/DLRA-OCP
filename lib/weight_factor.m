function [apply, unapply] = weight_factor(W)
%WEIGHT_FACTOR  Square-root factor of an SPD weight matrix.
%   [apply, unapply] = WEIGHT_FACTOR(W) returns two function handles
%       apply(X)   = C * X
%       unapply(X) = C \ X
%   where C'*C = W. W can be given as
%       * a vector of positive diagonal entries       (C = diag(sqrt(W)))
%       * a diagonal matrix, sparse or full (e.g. Ml) (C = diag(sqrt(diag(W))))
%       * a full SPD matrix (e.g. the consistent mass) (C = chol(W))
%   Used by WEIGHTED_QR and WEIGHTED_SVD to reduce the (M,Xi)-weighted
%   factorizations to Euclidean ones: X -> C_M X C_Xi'.

if isvector(W)
    w = full(W(:));
elseif isdiag(W)
    w = full(diag(W));
else
    C = chol(W);                     % W = C'*C, C upper triangular
    apply   = @(X) C * X;
    unapply = @(X) C \ X;
    return
end
if any(w <= 0)
    error('weight_factor:notSPD', 'Diagonal weights must be positive.');
end
s = sqrt(w);
apply   = @(X) s .* X;
unapply = @(X) X ./ s;
end