function X1 = explicit_rk(f, t0, t1, X0, tableau, nsub)
%EXPLICIT_RK  Explicit Runge-Kutta integration of a matrix ODE on [t0,t1].
%   X1 = EXPLICIT_RK(f, t0, t1, X0, tableau) integrates dX/dt = f(t,X) from
%   t0 to t1 starting from X0 with ONE step of the explicit Runge-Kutta
%   method given by TABLEAU (a name accepted by BUTCHER_TABLEAU, or a
%   struct with fields A, b, c). Default tableau: 'rk4'.
%   X1 = EXPLICIT_RK(f, t0, t1, X0, tableau, nsub) uses NSUB equal substeps.
%
%   X0 may be a matrix of any size (Nh x r, r x r, N x r, Nh x N, ...);
%   f(t,X) must return an array of the same size.
%
%   This is the workhorse behind MAKE_SUBSTEP_SOLVER for the K/L/S substeps
%   of the DLRA integrators (BUG, KSL, ...).

if nargin < 5 || isempty(tableau), tableau = 'rk4'; end
if nargin < 6 || isempty(nsub),    nsub = 1;        end
if ~isstruct(tableau), tableau = butcher_tableau(tableau); end

A = tableau.A; b = tableau.b; c = tableau.c;
s = numel(b);
h = (t1 - t0) / nsub;

X = X0;
t = t0;
for n = 1:nsub
    K = cell(s, 1);
    for i = 1:s
        Xi = X;
        for j = 1:i-1
            if A(i,j) ~= 0
                Xi = Xi + (h*A(i,j)) * K{j};
            end
        end
        K{i} = f(t + c(i)*h, Xi);
    end
    for i = 1:s
        if b(i) ~= 0
            X = X + (h*b(i)) * K{i};
        end
    end
    t = t + h;
end
X1 = X;
end