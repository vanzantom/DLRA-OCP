clear;clc
addpath('lib');addpath('solvers')
rng(0);
T=1;                                    %final time
N=100;                                  %matrix size
n1=50;      
I=eye(N,N);                        
Z=rand(N,N);                            %random matrix to estimate Lipschitz constant of projector

local=false;                             %true-> LTE, false->error at time t=T.


B1=rand(n1,N-n1);C1=rand(N-n1,N-n1);C1=0.5*(C1-C1');
W1=[zeros(n1,n1),-B1;B1',C1];
B2=rand(n1,N-n1);C2=rand(N-n1,N-n1);C2=0.5*(C2-C2');
W2=[zeros(n1,n1),-B2;B2',C2];
D=diag(2.^-(1:1:100));
Uf=@(t) expm(t*W1);
Vf=@(t) (expm(t*W2));
Sf=@(t) exp(t)*D;
Ufdot=@(t) W1*expm(t*W1);
Vfdot=@(t) W2*(expm(t*W2));
Sfdot=@(t) exp(t)*D;
Af=@(t) Uf(t)*Sf(t)*Vf(t)';
App=W1+eye(N,N);
Adot=@(t) App*Af(t)-Af(t)*W2;            %\dot A(t) = F(t,A(t))

rvec=[4;8;16;32];
Ntvec=2.^(4:1:9);
dtvec=T./Ntvec;
errfactorfact=zeros(length(rvec),length(dtvec));
errfactorproj=zeros(length(rvec),length(dtvec));
sigmavec=zeros(length(rvec),1);
cost=zeros(length(rvec),1);
effectivecurvature=zeros(length(rvec),Ntvec(end));

for i=1:length(rvec)
    r=rvec(i);

    
    Eperp = zeros(N,r); Eperp(r+1:2*r,:) = eye(r);
    E=zeros(N,r); E(1:r,:)=eye(r);
    U=@(t) Uf(t)*E; Udot=@(t) Ufdot(t)*E;
    V=@(t) Vf(t)*E; Vdot=@(t) Vfdot(t)*E;
    S=@(t) E'*Sf(t)*E; Sdot=@(t) E'*Sfdot(t)*E;
    Up=@(t) Uf(t)*Eperp; Vp=@(t) Vf(t)*Eperp;
    Updot=@(t) Ufdot(t)*Eperp; Vpdot=@(t) Vfdot(t)*Eperp;
    M=@(t) exp(t)*eye(r,r); Mdot=@(t) exp(t)*eye(r,r);
    A=@(t) U(t)*S(t)*V(t)';
    Adot=@(t) Udot(t)*S(t)*V(t)'+U(t)*Sdot(t)*V(t)'+U(t)*S(t)*Vdot(t)'; %A(t) is curve on \cM_r

    Up=@(t) Uf(t)*Eperp; 
    Vp=@(t) Vf(t)*Eperp;
    Updot=@(t) Ufdot(t)*Eperp; 
    Vpdot=@(t) Vfdot(t)*Eperp;
    

    func.U=U;func.Udot=Udot;func.V=V;func.Vdot=Vdot;func.Vp=Vp;
    func.Vpdot=Vpdot;func.Updot=Updot;func.Up=Up;
    func.M=M; func.Mdot=Mdot; func.S=S;
    for j=1:length(Ntvec)
        dt=T/Ntvec(j);
        if local==true
            Tfinal=dt;
        else
            Tfinal=T;
        end
        F=@(i) rhs_givendynamic(Adot);                                     
        Y=exact_dynamic(T,dt,func);                         
        
        Q=@(t) U(t)*M(t)*V(t)'+Up(t)*V(t)'+U(t)*Vp(t)';                 %tangent vector we want to approximate
        Q0.dU=Up(0);Q0.dS=M(0); Q0.dV=Vp(0);
 
        
        G=@(Q,t,Y) rhs_tangentODE(Q,t,Y,func);                          %projected Euler
        Qproj=solvecovariantODEProjection(@(Q,i)G(Q,(i-1)*dt,Y(i)),Q0,Y,Tfinal,dt,@(X,Z) projection(X,Z,eye(N),eye(N)));

        Gcont=@(Q,t) rhs_tangentODE(Q,t,Y,func);                        %integration of factors
        Qfact=integrate_factor_tangent_continuous(Gcont,Q0,Tfinal,dt,func);


        errfactorproj(i,j)=norm(tantodens(Qproj(end))-Q(Tfinal),'fro');
        errfactorfact(i,j)=norm(tantodens(Qfact(end))-Q(Tfinal),'fro');
        if isnan(errfactorproj(i,j))==1
            errfactorproj(i,j)=1e40;
        end
    end
    sigmavec(i)=min(svd(Y(end).S));
    X.U=U(0);X.V=V(0);X.S=S(0);                         %Estimate Lipschitz constant
    PXZ=projection(X,Z,I,I);
    Vptemp=Vp(0);
    X2.U=X.U;X2.S=X.S;X2.V=[X.V(:,1:r-1),Vptemp(:,1)];  %construction as in Kieri-Bart
    PYZ=projection(X2,Z,I,I);
    PXZ=X.U*PXZ.M*X.V'+PXZ.Up*X.V'+X.U*PXZ.Vp';
    PYZ=X2.U*PYZ.M*X2.V'+PYZ.Up*X2.V'+X2.U*PYZ.Vp';
    cost(i)=norm(PXZ-PYZ)/(norm(X2.U*X2.S*X2.V'-X.U*X.S*X.V')*norm(Z));

    Nt=Ntvec(end);
    for ll=1:Nt 
        effectivecurvature(i,ll)=norm(Udot((ll-1)*dt)*U((ll-1)*dt)'+U((ll-1)*dt)*Udot((ll-1)*dt)')+norm(Vdot((ll-1)*dt)*V((ll-1)*dt)'+V((ll-1)*dt)*Vdot((ll-1)*dt)');
    end

end
filename=sprintf('data_output/AnalysisCovariantODE.mat');
save(filename,'dtvec','errfactorproj','errfactorfact','sigmavec','cost','rvec','local','effectivecurvature');
