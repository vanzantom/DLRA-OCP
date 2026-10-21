function ConvergenceGradients(h,N,T)
%%CONVERGENCEGRADIENTS checks the convergence of the DLRA gradients with
%respect to \Delta t for several rank values
rng(0);
[geo,data,pars]=Setupfun(h,N,T,0.01);                                       %setup problem            
Ml=pars.Ml;clear pars;
Abasis=data.Abasis;chi=data.chi;Xi=data.Xi;Y0=data.Y0;forcing=data.forcing;

opts.substep='euler';                                                        %opts for time-integration
opts.n_boot=0;
F=@(z,i) affine_lowrank_rhs(Abasis,chi, Ml, Xi,forcing,z(:,i));
proj=@(Q,Y) projection(Y,Q,Ml,data.Xi);

dtvec=2.^-((7:1:12));
rvec=[4,8,16,32,64];

dtref=dtvec(end)/2;                                                         %compute Ref gradient on a finer time grid
data.dt=dtref; 
data.Nt=T/dtref;
[data]=generate_target(geo,data);                                           %compute target state (its size depends on \Delta t)
Yd=data.Yd;
tvec=linspace(0,1,data.Nt+1);tvec=tvec(1:end-1);                            %choice for the control
zinit=data.force*tvec;
S=@(z) solveODE(@(i)F(z,i),Y0,T,dtref,opts);
Fad=@(Y,i) affine_lowrank_rhs_adjoint(Abasis,chi, Ml, Xi,Y(:,:,i-1),Yd(:,i-1)); %i-1 is due to consistency with discretize-then-optimize
Sad=@(Y) solveODEback(@(i)Fad(Y,i),0*Y0,T,dtref,opts);
Jgrad=@(P,z) GradientFunctional(z,P,data);
Y=S(zinit);
P=Sad(Y);
gradRef=Jgrad(P,zinit);                                                     %grad Ref

errOtP=zeros(length(rvec),length(dtvec));
errPtO=zeros(length(rvec),length(dtvec));
errOtPAdjoint=zeros(length(rvec),length(dtvec));
errPtOAdjoint=zeros(length(rvec),length(dtvec));

for i=1:length(rvec)
    r=rvec(i);
    QT.dU=zeros(data.Nh,r);QT.dS=zeros(r,r);QT.dV=zeros(data.N,r);
    for j=1:length(dtvec)
        dt=dtvec(j);
        fprintf('Start computation for r=%d, dt=%f\n',r,dt);
        data.dt=dt; 
        data.Nt=T/dt;
        [data]=generate_target(geo,data);                                   %recompute target state.
        Yd=data.Yd;                                                         %redefine few functions that depend on \Delta t
        ErrMetric=@(z) Bochnorm(z,Ml,data);
        tvec=linspace(0,1,data.Nt+1);tvec=tvec(1:end-1);
        zinit=data.force*tvec;
        Jgrad=@(P,z) GradientFunctional(z,P,data); 
        
        SOtP=@(z) solveODEbug(@(i) F(z,i),Y0,T,dt,r,Ml,Xi,opts);            %gradient OtP approach
        Fadlr=@(Y,i) affine_lowrank_rhs_adjoint(Abasis,chi, Ml, Xi,Y(i),Yd(:,i));
        SadOtP=@(Y)  solveODEbugback(@(i)Fadlr(Y,i),0*Y0,T,dt,r,Ml,Xi);
        YOtP=SOtP(zinit);
        POtP=SadOtP(YOtP);
        gradOtP=Jgrad(POtP,zinit);   

        G=@(Y,z,i) rhs_covariant(Abasis, chi, Ml, Xi,forcing,dt,Y(i),Yd(:,i),z(:,i-1),i); %gradient PtO approach
        SPtO= @(z) solveODEbug(@(i)F(z,i),Y0,T,dt,r,Ml,Xi,opts);
        SadPtO=@(Y,z)  solvecovariantODEProjectionback(@(i) G(Y,z,i),QT,Y,T,dt,proj);
        PPtO=SadPtO(SPtO(zinit),zinit);
        gradPtO=Jgrad(PPtO,zinit);

        every=dt/dtref;
        errOtP(i,j)=sqrt(ErrMetric(gradRef(:,1:every:end)-gradOtP));
        errPtO(i,j)=sqrt(ErrMetric(gradRef(:,1:every:end)-gradPtO));
        errOtPAdjoint(i,j)=sqrt(ComputeNormAdjonts(P(:,:,1:every:end),POtP,Ml,Xi,data));
        errPtOAdjoint(i,j)=sqrt(ComputeNormAdjonts(P(:,:,1:every:end),PPtO,Ml,Xi,data));
    end
end
filename=sprintf('data_output/ConvergenceGradients_Nh_%d_N_%d%s.mat',data.Nh,data.N);
save(filename,'dtvec','rvec','errOtP','errPtO','errOtPAdjoint','errPtOAdjoint');