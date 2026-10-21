function A=exact_dynamic(T,dt,func)
%%EXACT_DYNAMIC returns an array of structures, containg the evaluation of
%%the time-dependent matrix function handles store in func.

Nt=T/dt;
U = func.U; V = func.V; S=func.S;

A = repmat(struct('U', [], 'S', [], 'V', []), Nt + 1, 1);
A(1).U=U(0); A(1).V=V(0); A(1).S=S(0);

for i=1:Nt
    A(i+1).U=U(i*dt); A(i+1).V=V(i*dt); A(i+1).S=S(i*dt);
end


