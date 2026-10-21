function [geo]=create_geometry_2D(h,basis_type)
%Geometry_Data receives the parameter h determining the mesh size and the
%basis_type specifing the FE space
%It generates the matrix geo which contain info about the mesh
    

    geo.regular=1;
    geo.h=h;
    geo.fun='@square';
    [P,T,Pb,Tb,E,Eb]=generate_mesh_2D(geo,basis_type);
    geo.P=P;
    geo.T=T;
    geo.Tb=Tb;
    geo.Pb=Pb;
    geo.E=E;
    geo.Eb=Eb;
    geo.Nv=size(geo.Pb,2); 
end

