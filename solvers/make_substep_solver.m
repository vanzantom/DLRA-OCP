function substep = make_substep_solver(method, varargin)
%MAKE_SUBSTEP_SOLVER  Build the integrator used for the K/L/S substeps.
%   substep = MAKE_SUBSTEP_SOLVER(method, ...) returns a function handle
%
%       X1 = substep(f, t0, t1, X0)
%
%   that integrates the matrix ODE dX/dt = f(t,X) from t0 to t1. This is
%   the single interface the DLRA integrators (BUG_STEP, ...) rely on, so
%   swapping the substep method never touches the integrator itself.
%
%   METHOD can be
%     * the name of an explicit Runge-Kutta method (see BUTCHER_TABLEAU):
%       'euler','heun','midpoint','rk3','rk4','rk38'.
%       Option 'nsub', k : split [t0,t1] into k equal RK substeps (default 1).
%     * the name of a MATLAB ODE solver: 'ode45','ode23','ode113','ode15s',...
%       Option 'options', odeset(...) (default RelTol=1e-8, AbsTol=1e-10).
%       The matrix is vectorised for the solver and reshaped back.
%     * a function handle @(f,t0,t1,X0) X1, returned unchanged. This is the
%       hook for user-written integrators (exponential Runge-Kutta, ...).
%
%   Examples
%       substep = make_substep_solver('rk4');
%       substep = make_substep_solver('heun', 'nsub', 4);
%       substep = make_substep_solver('ode45', 'options', odeset('RelTol',1e-6));
%       substep = make_substep_solver(@(f,t0,t1,X0) my_exp_rk(f,t0,t1,X0));

p = inputParser;
p.addParameter('nsub', 1);
p.addParameter('options', []);
p.parse(varargin{:});
nsub    = p.Results.nsub;
odeopts = p.Results.options;

if isa(method, 'function_handle')
    substep = method;
    return
end

switch lower(method)
    case {'euler','rk1','heun','rk2','midpoint','rk3','rk4','rk38'}
        tab = butcher_tableau(method);                 % built once
        substep = @(f, t0, t1, X0) explicit_rk(f, t0, t1, X0, tab, nsub);
    case {'ode45','ode23','ode113','ode15s','ode23s','ode23t','ode23tb','ode78','ode89'}
        solver = str2func(lower(method));
        if isempty(odeopts)
            odeopts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
        end
        substep = @(f, t0, t1, X0) matlab_ode_substep(solver, f, t0, t1, X0, odeopts);
    otherwise
        error('make_substep_solver:unknown', 'Unknown substep method ''%s''.', method);
end
end

function X1 = matlab_ode_substep(solver, f, t0, t1, X0, odeopts)
% Vectorise the matrix ODE for a MATLAB ODE suite solver and reshape back.
sz = size(X0);
vf = @(t, x) reshape(f(t, reshape(x, sz)), [], 1);
[~, xs] = solver(vf, [t0 t1], X0(:), odeopts);
X1 = reshape(xs(end, :).', sz);
end