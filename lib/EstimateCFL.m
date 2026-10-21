function CFL=EstimateCFL(data,M)
%%ESTIMATECFL estimates the CFL condition for the full-rank matrix ODE
%%vectorized

Abasis=data.Abasis;
m=data.Nh;
n=data.N;
mn=data.N*data.Nh;
Np=length(Abasis);
L=kron(speye(n),M\Abasis{1});
for k=1:Np-1
    D=diag(data.chi(:,k));
    L=L+kron(D,M\Abasis{k});
end
lambdamax=eigs(L,1,'largestabs');
CFL=2/lambdamax;
