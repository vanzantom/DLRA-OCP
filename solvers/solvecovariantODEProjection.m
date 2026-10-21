function Q=solvecovariantODEProjection(G,Q0,Y,T,dt,proj)
%SOLVEcovariantODEProjection integrates forward in time a covariant matrix ODE.
%   Q = solvecovariantODE_projection(F, Q0, Y, T, dt, proj) solves the covariant matrix ODE 
%           P_{Y(t)} \dot Q(t) = G(t,Q(t)), t\in (0,T]
%                Q(0) = Q0,
%   on the time grid t_i = i*dt, i = 0..Nt using a Projected Euler scheme.
%   The unknown Q(i) is an element of the tangent space at Y(i). It
%   consists of a structure containing the fields
%   Ubase (mxr), Vbase (nxr), dU (mxr), dS (rxr), dV (nxr) such that
%
%   Q(t_{i-1})=Q(i).Ubase*dS*Q(i).Vbase'+Q(i).dU*Q(i).Vbase'+Q(i).U*Q(i).dV'
%
%   Y(i) is also a structure such that Y(t_{i-1})=Y(i).U*Y(i).S*Y(i).V'.
%   
%   Since Q(t_i)\in T_{Y_{t_i)}M_r, it holds Y(i).U=Q(i).Ubase, Y(i).V=Q(i).Vbase 
%
%   G is a vector field s.t. G(Q,t,Y) returns an element of the tangent
%   space at Y(t).

                              
Nt=T/dt;
r=size(Y(1).S,1);

Q = repmat(struct('Ubase', [], 'Vbase', [],'dU', [], 'dS', [], 'dV', []), Nt+1, 1);      % du=Up, ds=M,dV=Vp.
Q(1).dU = Q0.dU; Q(1).dS = Q0.dS; Q(1).dV = Q0.dV;                                         
Q(1).Ubase=Y(1).U; Q(1).Vbase=Y(1).V;
for i = 1:Nt                                    % time stepping
    rhs=G(Q(i),i);                                %evaluate rhs in tangent format dU, dS, dV.
    Z.U=[Y(i).U,Q(i).dU+dt*rhs.dU];         %represent Q(i)+dt*F(Q(i)) in SVD format
    Z.V=[Y(i).V,Q(i).dV+dt*rhs.dV];
    Z.S=[Q(i).dS+dt*rhs.dS, eye(r); eye(r), zeros(r)];
    Zproj = proj(Y(i+1), Z);                    %project
    Q(i+1).dU = Zproj.Up; Q(i+1).dS = Zproj.M; Q(i+1).dV = Zproj.Vp; 
    Q(i+1).Ubase=Y(i+1).U; Q(i+1).Vbase=Y(i+1).V;
end