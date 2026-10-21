function [Q]=integrate_factor_tangent_continuous(F,Q0,T,dt,func)

Nt=T/dt;

U = func.U; Udot = func.Udot; 
V = func.V; Vdot = func.Vdot;


Q = repmat(struct('Ubase', [], 'Vbase', [],'dU', [], 'dS', [], 'dV', []), Nt+1, 1);      % du=Up, ds=M,dV=Vp.
Q(1).Ubase = U(0); Q(1).Vbase = U(0);
Q(1).dU = Q0.dU; Q(1).dS = Q0.dS; Q(1).dV = Q0.dV;

dt2 = dt / 2;
for i=1:Nt
      t=(i-1)*dt;

      %Euler
%     rhs=F(Q(i),told);
%     Q(i+1).dS=Q(i).dS +dt*(rhs.dS+Udot(t)'*Q(i).dU+Q(i).dV'*Vdot(t));
%     Q(i+1).dU=Q(i).dU +dt*(rhs.dU-Udot(t)*Q(i).dS-U(t)*(Udot(t)'*Q(i).dU));
%     Q(i+1).dV=Q(i).dV +dt*(rhs.dV-Vdot(t)*Q(i).dS-V(t)*(Vdot(t)'*Q(i).dV));
%     Q(i+1).Ubase = U(i*dt); 
%     Q(i+1).Vbase = V(i*dt);
  
    
    %RK4
    % --- Stage 1 (k1 at time t) ---
    k1 = eval_rhs(Q(i), t, F, U, V, Udot, Vdot);
    
    % --- Stage 2 (k2 at time t + dt/2) ---
    Q_stage2 = step_state(Q(i), k1, dt2, U, V, t + dt2);
    k2 = eval_rhs(Q_stage2, t + dt2, F, U, V, Udot, Vdot);
    
    % --- Stage 3 (k3 at time t + dt/2) ---
    Q_stage3 = step_state(Q(i), k2, dt2, U, V, t + dt2);
    k3 = eval_rhs(Q_stage3, t + dt2, F, U, V, Udot, Vdot);
    
    % --- Stage 4 (k4 at time t + dt) ---
    Q_stage4 = step_state(Q(i), k3, dt, U, V, t + dt);
    k4 = eval_rhs(Q_stage4, t + dt, F, U, V, Udot, Vdot);
    
    % --- Combine RK4 Stages ---
    Q(i+1).dS = Q(i).dS + (dt / 6) * (k1.dS + 2*k2.dS + 2*k3.dS + k4.dS);
    Q(i+1).dU = Q(i).dU + (dt / 6) * (k1.dU + 2*k2.dU + 2*k3.dU + k4.dU);
    Q(i+1).dV = Q(i).dV + (dt / 6) * (k1.dV + 2*k2.dV + 2*k3.dV + k4.dV);
    
    % Update base state for the next time index i*dt
    Q(i+1).Ubase = U(i*dt); 
    Q(i+1).Vbase = V(i*dt);
    tempU=(Q(i+1).Ubase'*Q(i+1).dU);
    tempV=(Q(i+1).Vbase'*Q(i+1).dV);
    Q(i+1).dU=Q(i+1).dU-Q(i+1).Ubase*tempU;
    Q(i+1).dV=Q(i+1).dV-Q(i+1).Vbase*tempV;
    Q(i+1).dS=Q(i+1).dS+tempU+tempV';
    

end
end


function rhs_tot = eval_rhs(Q_curr, t_curr, F, U, V, Udot, Vdot)
    % Call system RHS function
    rhs = F(Q_curr, t_curr);
    
    % Precompute time derivatives and base matrices at t_curr
    u  = U(t_curr);
    v  = V(t_curr);
    ud = Udot(t_curr);
    vd = Vdot(t_curr);
    
    % Calculate complete differential RHS
    rhs_tot.dS = rhs.dS + ud' * Q_curr.dU + Q_curr.dV' * vd;
    rhs_tot.dU = rhs.dU - ud * Q_curr.dS - u * (ud' * Q_curr.dU);
    rhs_tot.dV = rhs.dV - vd * Q_curr.dS - v * (vd' * Q_curr.dV);
end


%% Helper Function 2: Intermediate state update
function Q_out = step_state(Q_in, k_stage, step_size, U, V, t_next)
    Q_out.dS = Q_in.dS + step_size * k_stage.dS;
    Q_out.dU = Q_in.dU + step_size * k_stage.dU;
    Q_out.dV = Q_in.dV + step_size * k_stage.dV;
    
    % Assign corresponding base matrices at target stage time
    Q_out.Ubase = U(t_next);
    Q_out.Vbase = V(t_next);
end