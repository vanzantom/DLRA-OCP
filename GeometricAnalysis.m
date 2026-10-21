function []=GeometricAnalysis(h,N,T,dt)
%%EFFECTIVECURVATURE is a driver script to study a few quantities related to the DLRA approaches.
%It computes, for different ranks, the time evolution of
%   1) the effective curvature along the DLRA solution,
%   2) 1/sigma_r
%   3) the norm of \|\cP_Y^\perp F\|_{M,\Xi}
%   4) the norm of \|\cP_Q^\perp F_X\|_{M,\Xi}
%   5) the norm of \|\cP_Y^\perp F_X\|_{M,\Xi}
rng(0);
[~,data,pars]=Setupfun(h,N,T,dt);                 
Ml=pars.Ml;clear pars;
Abasis=data.Abasis;chi=data.chi;Xi=data.Xi;Yd=data.Yd;Y0=data.Y0;forcing=data.forcing;

F=@(z,i) affine_lowrank_rhs(Abasis,chi, Ml, Xi,forcing,z(:,i));
proj=@(Q,Y) projection(Y,Q,Ml,data.Xi);
tvec=linspace(0,1,data.Nt+1);tvec=tvec(1:end-1);                            %define a fixed control for every rank
zinit=data.force*tvec;


opts.n_boot=0;
rvec=[4,8,16,32];                                                           %rank values

kappaeff=zeros(data.Nt,length(rvec));                                       %Size of Y is data.Nt+1, but force field is evaluated Nt times
sigmar=zeros(data.Nt,length(rvec));
Fperp=zeros(data.Nt,length(rvec));
FxQperp=zeros(data.Nt,length(rvec));
FxYperp=zeros(data.Nt,length(rvec));
for i=1:length(rvec)
    r=rvec(i);
    fprintf('Start computation for r=%d\n',r);
    
  
    SPtO=@(z) solveODEbug(@(i) F(z,i),Y0,T,dt,r,Ml,Xi,opts);  
    Y=SPtO(zinit);                                                          %Compute DLRA state dynamics
    [k,sigma]=CurvatureEstimation(Y,@(i)F(zinit,i),dt,Ml,Xi,proj);        %Compute effective curvature, \sigma_r and norm of Fperp

    
    Fadlr=@(Y,i) affine_lowrank_rhs_adjoint(Abasis,chi, Ml, Xi,Y(i),Yd(:,i));
    SadOtP=@(Y)  solveODEbugback(@(i)Fadlr(Y,i),0*Y0,T,dt,r,Ml,Xi);
    POtP=SadOtP(Y);                                                         %Compute OtP adjoint
    
    
    G=@(Y,z,i) rhs_covariant(Abasis, chi, Ml, Xi,forcing,dt,Y(i),Yd(:,i),z(:,i-1),i);
    QT.dU=zeros(data.Nh,r);QT.dS=zeros(r,r);QT.dV=zeros(data.N,r);
    SadPtO=@(Y,z)  solvecovariantODEProjectionback(@(i) G(Y,z,i),QT,Y,T,dt,proj);
    PPtO=SadPtO(Y,zinit);                                                   %Compute PtO adjoint

    
    [Fp,FxQp,FxYp]=OrthogonalResiduals(Y,POtP,PPtO,@(i)F(zinit,i), @(i)Fadlr(Y,i),dt,Ml,Xi,proj);    %Compute orthogonal residuals
    
    kappaeff(:,i)=k;                                                        %save data
    sigmar(:,i)=sigma;
    Fperp(:,i)=Fp;
    FxQperp(:,i)=FxQp;
    FxYperp(:,i)=FxYp;
end
filename=sprintf('data_output/EffectiveCurvature_Nh_%d_N_%d_Deltat_%.2e.mat',data.Nh,data.N,dt);
save(filename,'tvec','rvec','kappaeff','sigmar','Fperp','FxQperp','FxYperp');
