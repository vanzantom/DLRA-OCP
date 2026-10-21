function [geo]=treat_boundary(geo)
%TREAT_BOUNDARY gets which nodes are on a Dirichlet boundary

nbn=size(geo.boundary,2);
v=ones(size(geo.Pb,2),1);
for k=1:nbn
    if geo.boundary(1,k)==-1% Dirichlet Data
        index=geo.boundary(3:end,k);% indexes of nodes which are on the edge boundary
        for j=1:length(index)
             v(index(j))=-1;
        end
    end
end
geo.int=find(v~=-1);%nodes not on the Dirichlet boundary
geo.out=find(v==-1); %nodes on the Dirichlet boundary

geo.Nh=length(geo.int);
end
