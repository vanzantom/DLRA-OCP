function [kappaeff,sigma]=CurvatureEstimation(Y,F,dt,M,Xi,proj)
%%CURVATUREESTIMATION computes the the time evolution of
%   1) the effective curvature along the DLRA solution,
%   2) 1/sigma_r
% Inputs: - DLRA solution Y
%         - Force Field F.
% 	      - Time step dt
% 	      - Matrices M and Xi defining the norm
%	      - Function handle proj to compute projections
% Outputs:
%         - kappaeff array of the time evolution of effective curvature
%	      - sigma array of time evolution of \sigma_r(Y).


sigma=zeros(length(Y)-1,1);
kappaeff=zeros(length(Y)-1,1);
LM=diag(sqrt(diag(M)));                    %Assumption: M and Xi are diagonal
LXi=diag(sqrt(diag(Xi)));
for i=1:length(Y)-1
    sigma(i)=min(svd(Y(i).S));
    rhs=F(i);
    Fval=rhs.Fstruct((i-1)*dt,Y(i));
    Ydot=proj(Fval,Y(i));       %Ydot=(M,Up,Vp): Ydot=UMV'+UpV'+UVp';
    DPuYdot=Ydot.Up*(Y(i).S\(Y(i).U'*M)) + Y(i).U*(Y(i).S'\(Ydot.Up'*M));
    DPvYdot=Ydot.Vp*(Y(i).S'\(Y(i).V'*Xi)) + Y(i).V*(Y(i).S\(Ydot.Vp'*Xi));
    kappaeff(i)=norm(LM*(DPuYdot/LM))+norm(LXi*(DPvYdot/LXi));
    
end
