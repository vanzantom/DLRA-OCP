%The script ``Setup" contains the parameters defining the model.
function [geo,data,pars]=Setupfun(h,N,T,dt)
    addpath('fem');
    addpath('lib');
    addpath('opt');
    addpath('data_input');
    addpath('solvers');
    rng(0);
                      % Problem definition
    data.h=h;      %mesh size
    data.T=T;         %Time interval
    data.dt=dt;   %time step
    data.Nt=data.T/data.dt;    
    data.p=4;         %# of random variables (p+1) basis matrices
    data.N=N;        %#samples   
    data.SC=0;        %SC=0 MC, SC=1, Stochastic collocation
    data.label=[-1,-1,-1,-1];       %Labels boundary conditions
    data.f=@(x,y) exp(x+y)*sin(2*pi*y); %Force term
    data.nu=10^(-3);  %penalization parameter;
    data.kcoeff=1*10^(-2);
    data.yd=@(x,y,t) t*sin(2*pi*x)*sin(2*pi*y);
    data.y0=@(x,y)   0*sin(2*pi*x)*cos(2*pi*x);
                      % Parameters for finite elements package  
    paraFEM.order_Gauss=2;  %Order Gauss quadrature.
    paraFEM.parallel=0;     %parallel assemble of stiffness matrices.
    paraFEM.basis_type=201; %P1
   
                      % Generate quadrature points and weights
    data=genquadpointsweights(data,data.p,data.N,data.SC);
    data.Xi=spdiags(data.weights,0,data.N,data.N);
                      % Mesh Generation
    filename=sprintf('data_input/geo_%d.mat',1/data.h);
    if(isfile(filename)~=1)
        GenerateGeometry(data.h,data,paraFEM);
    end
    load(filename,'geo');  
    data.Nh=length(geo.int); %size FEM space
                             % FEM Matrices Generation
    filename=sprintf('data_input/FEM_Matrices_%d.mat',1/data.h);
    %if(isfile(filename)~=1)
    generate_FEM_matrices(filename,geo,data,paraFEM); 
    %end
    [data,pars]=load_FEM_Matrices(filename,data);
    
    [data]=generate_target(geo,data);
    [pars.A]=GetMatrices(data);
    data.forcing=@(t) 0*exp(-t)*data.force;
    
   
    
                            
    