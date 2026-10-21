function output=tantodens(Q)
%%TANTODENS reiceves a tangent vector in a tangent format and returns its dense format
%
%Input: 
%       -Structure Q with fields Ubase,Vbase, dS, dU and dV
%        (Q=Ubase*dS*Vbase+ dUVbase' + UbasedV')
%Output:
%       -Tangent vector Q in its mxn dense format.

output=Q.Ubase*(Q.dS*Q.Vbase') + Q.dU*Q.Vbase' +  Q.Ubase*Q.dV';
end