function tab = butcher_tableau(name)
%BUTCHER_TABLEAU  Butcher tableau of a few explicit Runge-Kutta methods.
%   tab = BUTCHER_TABLEAU(name) returns a struct with fields
%       A : s x s strictly lower triangular RK matrix
%       b : s x 1 weights
%       c : s x 1 nodes (c = sum(A,2))
%   Available names:
%       'euler'    (= 'rk1')   explicit Euler, order 1
%       'heun'     (= 'rk2')   Heun,           order 2
%       'midpoint'             explicit midpoint, order 2
%       'rk3'                  Heun's 3-stage method, order 3
%       'rk4'                  classical Runge-Kutta, order 4
%       'rk38'                 3/8-rule, order 4
%   The same names are accepted by EXPLICIT_RK and MAKE_SUBSTEP_SOLVER.

switch lower(name)
    case {'euler','rk1'}
        A = 0;
        b = 1;
    case {'heun','rk2'}
        A = [0 0; 1 0];
        b = [1/2 1/2];
    case 'midpoint'
        A = [0 0; 1/2 0];
        b = [0 1];
    case 'rk3'
        A = [0 0 0; 1/3 0 0; 0 2/3 0];
        b = [1/4 0 3/4];
    case 'rk4'
        A = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];
        b = [1/6 1/3 1/3 1/6];
    case 'rk38'
        A = [0 0 0 0; 1/3 0 0 0; -1/3 1 0 0; 1 -1 1 0];
        b = [1/8 3/8 3/8 1/8];
    otherwise
        error('butcher_tableau:unknown', 'Unknown Runge-Kutta method ''%s''.', name);
end
tab.name = lower(name);
tab.A = A;
tab.b = b(:);
tab.c = sum(A, 2);
end