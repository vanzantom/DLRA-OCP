function [Fperp,FxQperp,FxYperp]=OrthogonalResiduals(Y,POtR,PRtO,F,Fx,dt,M,Xi,proj)
%%ORTHOGONALRESIDUALS computes the the time evolution of
%   1) the norm of \|\cP_Y^\perp F\|_{M,\Xi}
%   2) the norm of \|\cP_Q^\perp F_X\|_{M,\Xi}
%   3) the norm of \|\cP_Y^\perp F_X\|_{M,\Xi}
% Inputs: 
%        - DLRA solution Y
%	     - Force Field F.
% 	     - Time step dt
% 	     - Matrices M and Xi defining the norm
%	     - Function handle proj to compute projections
% Outputs:
%        - kappaeff array of the time evolution of effective curvature
%	     - sigma array for the time evolution of \sigma_r(Y).
%	     - Fperp array for the time evolution of \|\cP_Y^\perp F\|_{M,\Xi}
%        - FxQperp array for the time evolution of \|\cP_Q^\perp F_X\|_{M,\Xi}
%        - FxYperp array for the time evolution of \|\cP_Y^\perp F_X\|_{M,\Xi}


Fperp=zeros(length(Y)-1,1);
FxQperp=zeros(length(Y)-1,1);
FxYperp=zeros(length(Y)-1,1);

for i=1:length(Y)-1
    
    rhs=F(i);
    Fval=rhs.Fstruct((i-1)*dt,Y(i));
    Ydot=proj(Fval,Y(i)); 
    Fperp(i)=wnorm(Fval-Y(i).U*Ydot.M*Y(i).V'-Ydot.Up*Y(i).V'-Y(i).U*Ydot.Vp',M,Xi); %norm of \|\cP_Y^\perp F\|_{M,\Xi}
    
    rhsx=Fx(i);
    Fxval=rhsx.Fstruct((i-1)*dt,POtR(i));
    FxQtan=proj(Fxval,POtR(i));
    FxQperp(i)=wnorm(Fxval-POtR(i).U*FxQtan.M*POtR(i).V'-FxQtan.Up*POtR(i).V'-POtR(i).U*FxQtan.Vp',M,Xi); %norm of \|\cP_Q^\perp F_X\|_{M,\Xi}

    rhsx=Fx(i);
    Qdense=tantodens(PRtO(i));                                              %lazy implementation for the purpose of this analysis
    Fxval=rhsx.F((i-1)*dt,Qdense);
    FxYtan=proj(Fxval,Y(i));
    FxYperp(i)=wnorm(Fxval-Y(i).U*FxYtan.M*Y(i).V'-FxYtan.Up*Y(i).V'-Y(i).U*FxYtan.Vp',M,Xi); %norm of \|\cP_Y^\perp F_X\|_{M,\Xi}
   
end