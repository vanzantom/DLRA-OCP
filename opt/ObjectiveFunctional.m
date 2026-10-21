function J=ObjectiveFunctional(Y,z,Ms,data)
%%OBJECTIVEFUNCTIONAL evaluates the objective functional described in the
%%paper.
%
%Inputs: -control z
%        -state variable Y
%        -data structure
%Outputs:-scalar J, that is, the objective functiona value.

Nt=data.Nt;dt=data.dt;
J=0;

if isstruct(Y)              %If the state is passed in low-rank format I leverage structure to reduce C.C.
    e=ones(data.N,1);
    timeInt=0;
    for i=1:Nt              %loop over time steps
         Ydi=data.Yd(:,i);
         A=Y(i).U'*Ms*Y(i).U;
         B=Y(i).V'*data.Xi*Y(i).V;
         a=Y(i).U'*(Ms*Ydi);
         b=Y(i).V'*(data.Xi*e);
         T1=trace(B*(Y(i).S'*A*Y(i).S));
         T2=-2*(a'*Y(i).S*b);
         T3=(Ydi'*Ms*Ydi)*(e'*data.Xi*e);
         normsq=T1+T2+T3;
         timeInt=timeInt+normsq;
    end
    J=0.5*dt*timeInt;    
else
    avgsum=0;
    for j=1:data.N                                   %loop over samples                  
        timeInt=(Y(:,j,1)-data.Yd(:,1))'*Ms*(Y(:,j,1)-data.Yd(:,1));
        for i=2:Nt                                   %loop over time steps
             Yi=Y(:,j,i);
             Ydi=data.Yd(:,i);
             timeInt=timeInt+(Yi-Ydi)'*Ms*(Yi-Ydi);
        end
        avgsum=avgsum+0.5*dt*timeInt/data.N;
    end
    J=J+avgsum;
end

for i=1:Nt              %loop over time steps
    J=J+0.5*dt*data.nu*(z(:,i)'*Ms*z(:,i));    %penalization contribution
end