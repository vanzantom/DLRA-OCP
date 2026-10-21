function rhs=rhs_tangentODE(Q,t,Y,func)
%This function computes P_Y \dot Q assuming that Q(t)=U(t)M_0V(t).
% P_Y \dot Q = \dot U M_0 V^\top + U M_0 (\dot V)^\top. 

% RHS for projection Scheme for simple Q(t)=U(t)MV(t)'
% rhs.dS=zeros(size(Y.S));
% App=(Adot(t)*Q.Vbase*(Y.S\Q.dS));
% rhs.dU=App-Q.Ubase*(Q.Ubase'*App);
% App=(Adot(t)'*Q.Ubase*(Y.S'\Q.dS'));
% rhs.dV=App-Q.Vbase*(Q.Vbase'*App);

U = func.U; Udot = func.Udot; V = func.V;
Vdot = func.Vdot; Vpdot  = func.Vpdot;
Updot  = func.Updot;
M  = func.M;Mdot = func.Mdot; Up=func.Up; Vp=func.Vp;

rhs.dS=Mdot(t)+U(t)'*Updot(t)+Vpdot(t)'*V(t);
rhs.dU=Udot(t)*M(t)+Updot(t)-U(t)*(U(t)'*Updot(t));
rhs.dV=Vdot(t)*M(t)'+Vpdot(t)-V(t)*(V(t)'*Vpdot(t));

lambda=2;
rhs.dS=rhs.dS-lambda*(Q.dS-M(t));
rhs.dU=rhs.dU-lambda*(Q.dU-Up(t));
rhs.dV=rhs.dV-lambda*(Q.dV-Vp(t));

end
