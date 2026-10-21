function [U1, S1, V1, cache] = bug_step(U0, S0, V0, t0, t1, rhs, M, Xi, substep)
%BUG_STEP  One step of the fixed-rank BUG (basis update & Galerkin) integrator.
%   [U1,S1,V1] = BUG_STEP(U0,S0,V0, t0,t1, rhs, M, Xi, substep) advances the
%   rank-r approximation Y0 = U0*S0*V0' of dY/dt = F(t,Y) from t0 to t1 with
%   the "unconventional" integrator of Ceruti & Lubich (BIT 2022), written
%   in the (M,Xi)-weighted geometry of the paper:
%
%       U0'*M*U0 = I_r,   V0'*Xi*V0 = I_r   (and the same for U1, V1).
%
%   K-step   K(t0) = U0*S0,   dK/dt = F(t, K V0') Xi V0,       [U1,~] = QR_M (K(t1))
%   L-step   L(t0) = V0*S0',  dL/dt = F(t, U0 L')' M U0,       [V1,~] = QR_Xi(L(t1))
%   S-step   S(t0) = (U1' M U0) S0 (V0' Xi V1),
%            dS/dt = U1' M F(t, U1 S V1') Xi V1,               S1 = S(t1)
%
%   Inputs
%       rhs      struct with factories rhs.K(V), rhs.L(U), rhs.S(U,V), each
%                returning a handle @(t,X) (see AFFINE_LOWRANK_RHS)
%       M, Xi    SPD weights (vector, diagonal or full matrix)
%       substep  handle X1 = substep(f,t0,t1,X0) (see MAKE_SUBSTEP_SOLVER)
%   Output
%       cache    struct with the QR factors RK, RL, the reduced blocks
%                Mmat = U1'MU0, Nmat = V1'XiV0 and S_tilde, which are what
%                the discrete adjoint (reverse mode) of the step needs.

% K-step: update the left basis, old right basis frozen
K1 = substep(rhs.K(V0), t0, t1, U0 * S0);
[U1, RK] = weighted_qr(K1, M);

% L-step: update the right basis, old left basis frozen
L1 = substep(rhs.L(U0), t0, t1, V0 * S0');
[V1, RL] = weighted_qr(L1, Xi);

% S-step: Galerkin step on the new bases
Mmat = U1' * wmul(M, U0);            % U1' M  U0   (r x r)
Nmat = V1' * wmul(Xi, V0);           % V1' Xi V0   (r x r)
S_tilde = Mmat * S0 * Nmat';
S1 = substep(rhs.S(U1, V1), t0, t1, S_tilde);

if nargout > 3
    cache = struct('RK', RK, 'RL', RL, 'Mmat', Mmat, 'Nmat', Nmat, 'S_tilde', S_tilde);
end
end

function Y = wmul(W, X)
% W*X for W given as a vector of diagonal entries or as a matrix
if isvector(W)
    Y = W(:) .* X;
else
    Y = W * X;
end
end