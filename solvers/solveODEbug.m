function [Y,info] = solveODEbug(F,Y0,T,dt,r,M,Xi,opts)
%SOLVEODEBUG integrates forward in time an ODE with 
%the fixed-rank r BUG integrator.
%   Y = solveODEbug(F,Y0,T,dt,r,opts) is the low-rank counterpart of solveODE to solve
%           \dot Y(t) = F(t,Y(t)), t\in (0,T]
%                Y(0) = Y0,
%   on the time grid t_i = i*dt, i = 0..Nt, 
%   The state is kept in factorized form
%
%       Y(i).U * Y(i).S * Y(i).V'  ~  Y(t_{i-1}),    i = 1, ..., Nt+1,
%
%   with Y(i).U (Nh x r), Y(i).S (r x r), Y(i).V (N x r) and the gauge
%   U'*M*U = I_r, V'*Xi*V = I_r, Xi = diag(data.weights). The dense state
%   at t_{i-1} is recovered as Y(i).U*Y(i).S*Y(i).V' (Nh x N, cost O(Nh N r)).
%
%   The vector field is specificied by the function handle F which can be
%   evaluated for every i =0..Nt.
%
%   [Y, info] = solveODEbug(...) also returns info.n_boot and the
%   per-step caches info.cache(i) of BUG_STEP (QR factors etc.) for i >
%   n_boot
%
%   opts (all optional):
%     .substep  'euler' (default) | 'heun' | 'midpoint' | 'rk3' | 'rk4' |
%               'rk38' | 'ode45' | ... | function handle @(f,t0,t1,X0),
%               see MAKE_SUBSTEP_SOLVER. Also used for the bootstrap steps.
%     .nsub     number of RK substeps inside each K/L/S step (default 1)
%     .odeopts  odeset options when .substep is a MATLAB ODE solver
%     .n_boot   number of bootstrap steps (default min(r+2, 10), see below)
%
%   Bootstrap. The initial condition Y0 = 0 (and more generally a rank-
%   deficient Y0) does not belong to the rank-r manifold. To address this,
%   the first n_boot steps are taken with the DENSE scheme
%   (same substep method applied to rhs.F) and each dense state is
%   truncated to rank r by WEIGHTED_SVD; the BUG chain then starts from the
%   truncation of the n_boot-th state. Set opts.n_boot = 0 to start BUG
%   directly from the rank-r truncation of Y0.

if nargin < 8, opts = struct(); end
getopt = @(name, default) getfield_or(opts, name, default);
Nt=T/dt; Nh = size(Y0,1); N=size(Y0,2);

if r > min(Nh, N)
    error('solve_state_bug:rank', 'rank r=%d must be <= min(Nh,N)=%d', r, min(Nh, N));
end

n_boot  = min(getopt('n_boot', min(r + 2, 10)), Nt);

substep = make_substep_solver(getopt('substep', 'euler'), ...
                              'nsub', getopt('nsub', 1), ...
                              'options', getopt('odeopts', []));

Y = repmat(struct('U', [], 'S', [], 'V', []), Nt + 1, 1);

                                            % bootstrap: dense steps, each truncated to rank r                                               
[Y(1).U, Y(1).S, Y(1).V] = weighted_svd(Y0, r, M, Xi);
for i = 1:n_boot
    rhs = F(i);
    Y0  = substep(rhs.F, (i-1)*dt, i*dt, Y0);
    [Y(i+1).U, Y(i+1).S, Y(i+1).V] = weighted_svd(Y0, r, M, Xi);
end

% ---- BUG steps ------------------------------------------------------------
U = Y(n_boot+1).U;  S = Y(n_boot+1).S;  V = Y(n_boot+1).V;
if nargout > 1
    info.n_boot = n_boot;
    info.cache  = repmat(struct('RK',[],'RL',[],'Mmat',[],'Nmat',[],'S_tilde',[]), Nt, 1);
end
for i = n_boot+1:Nt
    rhs = F(i);
    if nargout > 1
        [U, S, V, info.cache(i)] = bug_step(U, S, V, (i-1)*dt, i*dt, rhs, M, Xi, substep);
    else
        [U, S, V] = bug_step(U, S, V, (i-1)*dt, i*dt, rhs, M, Xi, substep);
    end
    Y(i+1).U = U;  Y(i+1).S = S;  Y(i+1).V = V;
end

function v = getfield_or(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end
end