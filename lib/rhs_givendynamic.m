function rhs = rhs_givendynamic(Adot)
%rhs_givendynamic compute the right-hand side (dense and K/L/S) for the matrix
%ODE    \dot A(t) = F(t,A).
% The code receives as input the time derivative of A, Adot, thus F(t,A).
%   The returned struct contains four function handles:
%       rhs.F(t,Y)   dense right-hand side, Y is Nh x N
%       rhs.K(V)     handle @(t,K) F(t, K V') Xi V          (K-step, K is Nh x r)
%       rhs.L(U)     handle @(t,L) F(t, U L')' M U          (L-step, L is N x r)
%       rhs.S(U,V)   handle @(t,S) U' M F(t, U S V') Xi V   (S-step, S is r x r)

rhs.F = @(t, Y) Adot(t);
rhs.K = @(V)    makeK(V,Adot);
rhs.L = @(U)    makeL(U,Adot) ;
rhs.S = @(U, V) makeS(U,V,Adot);
end

function f=makeK(V,Adot)
f=@(t,K) Adot(t)*V;
end
function f=makeL(U,Adot)
f=@(t,L) Adot(t)'*U;
end
function f=makeS(U,V,Adot)
f=@(t,S) U'*Adot(t)*V;
end
