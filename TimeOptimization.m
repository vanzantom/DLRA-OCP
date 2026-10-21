function TimeOptimization(h,N,T,dt)
%TIMEOPTIMIZATION solves the minimization problem with the full-rank and the
%two DLRA approaches.
%The goal is to collect data to compare the computational times.
rng(0);
[geo,data,pars]=Setupfun(h,N,T,dt);                 
Ml=pars.Ml;clear pars;
Abasis=data.Abasis;chi=data.chi;Xi=data.Xi;Yd=data.Yd;Y0=data.Y0;forcing=data.forcing;

maxit=30;
opts.substep='euler';
opts.n_boot=0;
F=@(z,i) affine_lowrank_rhs(Abasis,chi, Ml, Xi,forcing,z(:,i));
proj=@(Q,Y) projection(Y,Q,Ml,data.Xi);
ErrMetric=@(z) Bochnorm(z,Ml,data);
zinit=zeros(geo.Nh,data.Nt);

S=@(z) solveODE(@(i)F(z,i),Y0,T,data.dt,opts);
Fad=@(Y,i) affine_lowrank_rhs_adjoint(Abasis,chi, Ml, Xi,Y(:,:,i-1),Yd(:,i-1));  %i-1 is due to consistency with discretize-then-optimize
Sad=@(Y,z) solveODEback(@(i)Fad(Y,i),0*Y0,T,dt,opts);
Jgrad=@(P,z) GradientFunctional(z,P,data);
J=@(Y,z) ObjectiveFunctional(Y,z,Ml,data);

tic;
[z,iter,errvec,Jvec]=myGD(zinit,J,Jgrad,S,Sad,ErrMetric,maxit,1e-6,1/(4*data.nu));
timefull=toc;

rvec=[4,8,16,32];

timeOtP=zeros(length(rvec),1);
timePtO=zeros(length(rvec),1);
JvecOtP=zeros(length(rvec),maxit+1);
JvecPtO=zeros(length(rvec),maxit+1);
zPtO=cell(length(rvec),1);
zOtP=cell(length(rvec),1);
errvecOtP=zeros(length(rvec),maxit);
errvecPtO=zeros(length(rvec),maxit);


for i=1:length(rvec)
    r=rvec(i);
    fprintf('Start computation for r=%d, dt=%f\n',r);
    QT.dU=zeros(data.Nh,r);QT.dS=zeros(r,r);QT.dV=zeros(data.N,r);
    
    SOtP=@(z) solveODEbug(@(i) F(z,i),Y0,T,dt,r,Ml,Xi,opts);      
    Fadlr=@(Y,i) affine_lowrank_rhs_adjoint(Abasis,chi, Ml, Xi,Y(i),Yd(:,i));
    SadOtP=@(Y,z)  solveODEbugback(@(i)Fadlr(Y,i),0*Y0,T,dt,r,Ml,Xi);
    
    tic;
    [zOtPtemp,iterOtP,errvecOtPtemp,JvecOtPtemp]=myGD(zinit,J,Jgrad,SOtP,SadOtP,ErrMetric,maxit,1e-6,1/(4*data.nu));
	timeOtP(i)=toc;
    JvecOtP(i,:)=JvecOtPtemp';
    zOtP{i}=zOtPtemp;
    errvecOtP(i,:)=errvecOtPtemp';
    
    G=@(Y,z,i) rhs_covariant(Abasis, chi, Ml, Xi,forcing,dt,Y(i),Yd(:,i),z(:,i-1),i);
    SPtO= @(z) solveODEbug(@(i)F(z,i),Y0,T,dt,r,Ml,Xi,opts);
    SadPtO=@(Y,z)  solvecovariantODEProjectionback(@(i) G(Y,z,i),QT,Y,T,dt,proj);
    
    tic;
    [zPtOtemp,iterPtO,errvecPtOtemp,JvecPtOtemp]=myGD(zinit,J,Jgrad,SPtO,SadPtO,ErrMetric,maxit,1e-6,1/(4*data.nu));
	timePtO(i)=toc;
    JvecPtO(i,:)=JvecPtOtemp';
    zPtO{i}=zPtOtemp;
    errvecPtO(i,:)=errvecPtOtemp';
end

Full.iter=iter; Full.errvec=errvec; Full.Jvec=Jvec; Full.time=timefull;
OtP.iter=iterOtP; OtP.err=errvecOtP; OtP.Jvec=JvecOtP; OtP.time=timeOtP;
PtO.iter=iterPtO; PtO.err=errvecPtO; PtO.Jvec=JvecPtO; PtO.time=timePtO;
%Full.z=z; OtP.z=zOtP; PtO.z=zPtO;      %avoid storing controls if not needed

filename=sprintf('data_output/TimeOptimization_Nh_%d_N_%d%s.mat',data.Nh,data.N);
save(filename,'rvec','OtP','PtO','Full');