function Y=solveODEback(F,Y_T,T,dt,opts)
%SOLVEODE integrates backward in time a matrix ODE.
%   Y = solveODEback(F, Y0, T, dt, opts) solves the matrix ODE 
%           \dot Y(t) = F(t,Y(t)), t\in (0,T]
%                Y(T) = Y_T,
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
Y = zeros(size(Y_T,1),size(Y_T,2),Nt+1);
Y(:,:,Nt+1) = Y_T;                                              
for i = Nt:-1:1        % time stepping
    rhs = F(i+1);          
    Y(:,:,i)  = substep(rhs.F, (i-1)*dt, i*dt, Y(:,:,i+1));
end

end 

function v = getfield_or(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end