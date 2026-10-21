function Y=solveODE(F,Y0,T,dt,opts)
%SOLVEODE integrates forward in time a matrix ODE.
%   Y = solveODE(F, Y0, T, dt, opts) solves the matrix ODE 
%           \dot Y(t) = F(t,Y(t)), t\in (0,T]
%                Y(0) = Y0,
%   on the time grid t_i = i*dt, i = 0..Nt, 
%   The matrix Y(:,:,i) is of size Nh x N and approximates Y(t) at time t_{i-1}.
%   The vector field is specificied by the function handle F which can be
%   evaluated for every i =0..Nt.

%   opts (all optional):
%     .substep  'euler' (default) | 'heun' | 'midpoint' | 'rk3' | 'rk4' |
%               'rk38' | 'ode45' | ... | function handle @(f,t0,t1,X0),
%               see MAKE_SUBSTEP_SOLVER.
%     .odeopts  odeset options when .substep is a MATLAB ODE solver

if nargin < 5, opts = struct(); end
getopt = @(name, default) getfield_or(opts, name, default);

substep = make_substep_solver(getopt('substep', 'euler'), ...
                              'options', getopt('odeopts', []));
                              
Nt=T/dt;
Y = zeros(size(Y0,1),size(Y0,2),Nt+1);
Y(:,:,1) = Y0;                                              
for i = 1:Nt        % time stepping
    rhs = F(i);
    Y(:,:,i+1)  = substep(rhs.F, (i-1)*dt, i*dt, Y(:,:,i));
end

end

function v = getfield_or(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end