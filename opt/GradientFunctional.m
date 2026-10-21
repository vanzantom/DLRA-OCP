function grad=GradientFunctional(z,P,data)
%%GRADIENTFUNCTIONAL evaluates the gradient: \nablaJ(z)=nu*z-B\star P.
%
%Inputs: -control z
%        -adjoint variable P
%        -data structure
%Outputs:-grad, a matrix Nh\times Nt, Nh size FE space, Nt number of time discretization points

Nt=data.Nt;
grad=zeros(data.Nh,data.Nt);                

if isstruct(P)                                           %distinguish whether P is
    if numel(fieldnames(P))==3                           %a point on the manifold
        for i=1:Nt
            w=mean(P(i+1).V,1)';   % r x 1
            Pmean=P(i+1).U * (P(i+1).S * w);
            grad(:,i)=data.nu*z(:,i)-Pmean;
        end
    elseif numel(fieldnames(P))==5                       %an element of the tangent space
       for i=1:Nt
            w=mean(P(i+1).Vbase,1)';   % r x 1
            Pmean=P(i+1).Ubase * (P(i+1).dS * w);
            w=mean(P(i+1).Vbase,1)';
            Pmean=Pmean+P(i+1).dU * w;
            w=mean(P(i+1).dV,1)';
            Pmean=Pmean+P(i+1).Ubase * w;
            
            grad(:,i)=data.nu*z(:,i)-Pmean;
        end
        
    else
        disp('Error in the format of P!');
        exit;
    end
else                                                     %a matrix
    for i=1:Nt                                          %loop over time steps
        grad(:,i)=data.nu*z(:,i)-mean(P(:,:,i+1),2);    %vectorized mean
    end
end