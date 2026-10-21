function [data]=generate_target(geo,data)


Yd=zeros(data.Nh,data.Nt+1);                  % target state
for i=1:data.Nt+1                             %Yd{i} is a matrix Ydm\in R^{Nh x Nt}
    ytemp=zeros(data.Nh,1);                   %assumption yd is deterministic (can be removed)
    for k=1:length(geo.int)
        ytemp(k)=data.yd(geo.Pb(1,geo.int(k)),geo.Pb(2,geo.int(k)),i*data.dt);
    end
    Yd(:,i)=ytemp;
end
data.Yd=Yd;

