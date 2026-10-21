function rhs = affine_lowrank_rhs_adjoint(Abasis, chi, M, Xi,Y,Yd)
%AFFINE_LOWRANK_RHS_adjoint  compute the right-hand side (dense and K/L/S) of the adjoint equation.
%   rhs = AFFINE_LOWRANK_RHS_adjoint(Abasis, chi, M, Xi, g) computes the vector field
%       for the the matrix ODE
%
%       dY/dt=F(t,Y)=-M^{-1}(sum_{k=0}^{N_p} A_k Y D_k ) + G(t) 
%
%       where
%       Abasis : cell {A_0, A_1, ..., A_{N_p}} of Nh x Nh (sparse) stiffness matrices
%       chi    : N x p matrix of samples, D_k = diag(chi(:,k))
%       M      : Nh x Nh (possibly lumped) mass matrix
%       Xi     : N x N diagonal matrix containing quadrature weights
%       G      : function returning force term G(t) of size Nh x N in its low rank format
%                G(t).Phi*G(t).Psi, G(t)= ones(1,N)*Yd(t)-Y(t) (both Y, Y_d constant
%                on the time subinterval)
%   Column j of F is -M^{-1} A(omega_j) y_j + G(:,j), with A(omega_j) =
%   A_0 + sum_k chi(j,k) A_k, G=g.Phi*g.Psi'
%
%   The returned struct contains four function handles:
%       rhs.F(t,Y)   dense right-hand side, Y is Nh x N
%       rhs.K(V)     handle @(t,K) F(t, K V') Xi V          (K-step, K is Nh x r)
%       rhs.L(U)     handle @(t,L) F(t, U L')' M U          (L-step, L is N x r)
%       rhs.S(U,V)   handle @(t,S) U' M F(t, U S V') Xi V   (S-step, S is r x r)
%   with the gauge U'*M*U = I_r, V'*Xi*V = I_r. The reduced matrices
%       B_k = V' D_k Xi V  (r x r),   C_k = U' A_k U  (r x r),
%       v   = V' Xi 1_N     (r x 1),   MU  = M U        (Nh x r)
%   are precomputed once per factory call, so each evaluation of a substep
%   right-hand side costs O((N_p+1)(Nh + N) r^2) plus one sparse M-solve
%   (K-step only), and never forms an Nh x N matrix.

if isvector(Xi), w = full(Xi(:)); else, w = full(diag(Xi)); end
if nargin < 6 Yd=[]; end
if nargin < 5  Y=[]; end         
N=length(w);

if isstruct(Y)==0
    [U,S,V]=svd(Y);
    clear Y;
    Y.U=U; Y.S=S; Y.V=V;
end
Ylr.Phi=Y.U*Y.S;
Ylr.Psi=Y.V;
g.Phi = @(t) [Yd, -Ylr.Phi];
g.Psi=@(t) [ones(N,1),Ylr.Psi];           


rhs.F = @(t, Y) -(M \ applyA(Y, Abasis, chi)) + g.Phi(t)*g.Psi(t)';
rhs.Fstruct = @(t, Y) -(M \ applyA_struct(Y, Abasis, chi)) + g.Phi(t)*g.Psi(t)';
rhs.K = @(V)    makeK(V, Abasis, chi, M, w, g);
rhs.L = @(U)    makeL(U, Abasis, chi, M, g);
rhs.S = @(U, V) makeS(U, V, Abasis, chi, M, w, g);
end

function AY = applyA(Y, Abasis, chi)               % dense
% A_0 Y + sum_k A_k Y D_k, D_k diagonal
AY = Abasis{1} * Y;
for k = 1:size(chi, 2)
    AY = AY + (Abasis{k+1} * Y) .* chi(:, k)';
end
end

function AY = applyA_struct(Y, Abasis, chi)        % input low-rank format. inefficient implementation (used for GeometricAnalyses.m)
% A_0 Y + sum_k A_k Y D_k, D_k diagonal
AY = Abasis{1} * (Y.U*Y.S*Y.V');
for k = 1:size(chi, 2)
    AY = AY + (Abasis{k+1} * (Y.U*Y.S*Y.V')) .* chi(:, k)';
end
end

% ------------------------------------------------------------------ K-step

function f = makeK(V, Abasis, chi, M, w, g) % K step
[B, v] = reducedB(V, chi, w);
f = @(t, K) -(M \ sumAKB(K, B, Abasis)) + g.Phi(t) * (g.Psi(t)'*v);
end

function [B, v] = reducedB(V, chi, w)
% B_k = V' D_k Xi V (B_0 = V' Xi V = I_r by the gauge)
XiV = w .* V;
B = cell(size(chi, 2) + 1, 1);
B{1} = V' * XiV;
for k = 1:size(chi, 2)
    B{k+1} = (chi(:, k) .* V)' * XiV;
end
v = XiV;
end

function out = sumAKB(K, B, Abasis)
out = (Abasis{1} * K) * B{1};
for k = 2:numel(Abasis)
    out = out + (Abasis{k} * K) * B{k};
end
end

function f = makeL(U, Abasis, chi, M, g)
C  = reducedC(U, Abasis);
MU = M * U;
f = @(t, L) -sumLC(L, C, chi) + g.Psi(t) * (MU' * g.Phi(t))';
end

function C = reducedC(U, Abasis)
% C_k = U' A_k U (= U' M (M^{-1} A_k) U, the mass inverse cancels)
C = cell(numel(Abasis), 1);
for k = 1:numel(Abasis)
    C{k} = U' * (Abasis{k} * U);
end
end

function out = sumLC(L, C, chi)
out = L * C{1}';
for k = 1:size(chi, 2)
    out = out + chi(:, k) .* (L * C{k+1}');      % D_k L C_k'
end
end

% ------------------------------------------------------------------ S-step
function f = makeS(U, V, Abasis, chi, M, w, g)
C  = reducedC(U, Abasis);
MU = M * U;
[B, v] = reducedB(V, chi, w);
f = @(t, S) -sumCSB(S, C, B) + (MU' * g.Phi(t)) * (g.Psi(t)'*v);
end

function out = sumCSB(S, C, B)
out = C{1} * S * B{1};
for k = 2:numel(C)
    out = out + C{k} * S * B{k};
end
end