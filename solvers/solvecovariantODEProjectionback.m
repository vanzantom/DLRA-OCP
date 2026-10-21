function Q=solvecovariantODEProjectionback(G,QT,Y,T,dt,proj)
%SOLVEcovariantODEProjectionback integrates backward in time a covariant matrix ODE.
%   Q = solvecovariantODE_projection(F, QT, Y, T, dt, proj) solves the covariant matrix ODE 
%           -P_{Y(t)} \dot Q(t) = G(t,Q(t)), t\in [0,T)
%                Q(T) = QT,
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
%   G is a vector field s.t. G(Q,t,Y) returns an element of the tangent
%   space at Y(t).

                              
Nt=T/dt;
r=size(Y(1).S,1);

Q = repmat(struct('Ubase','Vbase','dU', [], 'dS', [], 'dV', []), Nt+1, 1);      % du=Up, ds=M,dV=Vp.
Q(Nt+1).dU = QT.dU; Q(Nt+1).dS = QT.dS; Q(Nt+1).dV = QT.dV;                                         
Q(Nt+1).Ubase=Y(Nt+1).U; Q(Nt+1).Vbase=Y(Nt+1).V;

for i = Nt:-1:1                                    % time stepping
    rhs=G(i+1);                                    %evaluate rhs in tangent format dU, dS, dV.
    rhs_tan=rhs(Q(i+1));
    Z.U=[Y(i+1).U,Q(i+1).dU+dt*rhs_tan.dU];       %represent Q(i+1)+dt*G(t_{i+1},Q(i+1)) in SVD format
    Z.V=[Y(i+1).V,Q(i+1).dV+dt*rhs_tan.dV];
    Z.S=[Q(i+1).dS+dt*rhs_tan.dS, eye(r); eye(r), zeros(r)];
    Zproj = proj(Z, Y(i));                          %project
    Q(i).dU = Zproj.Up; Q(i).dS = Zproj.M; Q(i).dV = Zproj.Vp; 
    Q(i).Ubase=Y(i).U; Q(i).Vbase=Y(i).V;
end