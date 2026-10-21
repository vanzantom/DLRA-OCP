function  [z,iter,errvec,Jvec]=myGD(z,J,Jgrad,S,Sad,Bochnorm,maxit,tol,alphainit)
%%MYGD implements a gradient descent method with Armijo backtracking
%
%Inputs:
%       -z initial control
%       -J objective functional
%       -J functional gradient
%       -S solution operator
%       -Sad adjoint operator
%       -Bochnorm function to compute norms of errors
%       -maxit maximum number of iteration
%       -tol tolerance on the gradient norm
%       -alphainit initial step size length
%Outputs:
%       -z final control
%       -iter array of size equal to maximum number of iterations
%       -errvec array containing gradient norm at each iteration
%       -Jvec array containing objective functional values at each it.

rho=0.5;
gamma=1e-4;
errvec=zeros(maxit,1);
Jvec=zeros(maxit+1,1); 
Y=S(z);
Jval=J(Y,z);
Jvec(1)=Jval;
k=1;
flag=0;                           %flag will be 1 if tolerance reached
while(k<=maxit && flag==0)
    alpha=alphainit;               
    P=Sad(Y,z);
    grad=Jgrad(P,z);
    znew=z-alpha*grad;
    Ynew=S(znew);
    Jvalnew=J(Ynew,znew);
    j=1;
    flagls=0;
    while(Jvalnew-Jval>-gamma/alpha*Bochnorm(z-znew)  && j<20) %descent condition
        alpha=alpha*rho^j;
        znew=z-alpha*grad;
        Ynew=S(znew);
        Jvalnew=J(Ynew,znew);
        j=j+1;
        flagls=1;
    end
    err=sqrt(Bochnorm(grad));         %norm in L^2([0,T],L^2(D))
    if err<tol
        flag=1;
    end
    errvec(k)=err;
    z=znew;
    Jvec(k+1)=Jvalnew;
    Jval = Jvalnew;
    Y=Ynew;
    fprintf('It: %d|| J=%.7f || norm Grad =%.7f || # L.S.: %d\n',k,Jvalnew,err,j-1);
    k=k+1;
    if(j==20)
        fprintf('Line search did not converge...exiting\n');
        break;
    end
end
errvec=errvec(1:k-1);
Jvec=Jvec(1:k); 
iter=k-1;