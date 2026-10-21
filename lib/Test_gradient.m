function []=Test_gradient(data,S,Sad,Jgrad,J,Ms,z)
%Check analytic gradient with finite difference approximation

if nargin<7
    z=rand(data.Nh,data.Nt);
end
dz=rand(data.Nh,data.Nt);                     %random perturbation
Y=S(z);
P=Sad(Y,z);
grad=Jgrad(P,z);
Jval=J(Y,z);

                                              %First test: finite differences
% epsi=1;
% for j=1:10       
%     epsivec(j)=epsi;
%     Yplus=S(z+epsi*dz);
%     Jvalplus=J(Yplus,z+epsi*dz);
%     test(j) = abs( ActionGradientScalarProduct(grad,Ms,data.T,dz) - (Jvalplus-Jval)/(epsi));
%     epsi=epsi/10;
%     [epsi test(j)]
%     loglog(epsivec,test(j))
% end
% figure(1)
% loglog(epsivec,test,'b-o');hold on
% loglog(epsivec,epsivec,'k-o');
                                              %second test: Taylor
epsi=1;
for j=1:10      
    epsivec(j)=epsi;
    Yplus=S(z+epsi*dz);
    Jvalplus=J(Yplus,z+epsi*dz);
    test2(j) = abs(Jvalplus-Jval-epsi*ActionGradientScalarProduct(grad,Ms,data.T,dz));
    epsi=epsi/10;
    [epsi test2(j)];
end
figure(2)
loglog(epsivec,test2,'b-o');hold on
loglog(epsivec,epsivec.^2,'k-o');
end


function res=ActionGradientScalarProduct(grad,Ms,T,dz)

Nt=size(grad,2);
res=0;
dt=T/Nt;
for i=1:Nt
    res=res+dt*(grad(:,i)'*Ms*dz(:,i));
end

end