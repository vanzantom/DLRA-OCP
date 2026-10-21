%-------------------------------------------------------------------------
% generate_mesh_2D generates the mesh on a geometry specified by the geometric function geo.fun with
% maximum mesh size geo.h.
% generate_mesh_2D receives
% geo: structure containing a function handle describing the geometry and
% the meshsize geo.h and a flag for a regular/non regular mesh.
% basis_type: type of FEM space (linear,quadratic etc..).
% generate_mesh_2D returns:
% P: vertex matrix containing position of vertices of the mesh. Eg: P(:,1)=[X_1;Y_1]
% T: element matrix containing indeces vertices of the elements plus subdomain to which element k belongs. 
%     Eg: T(:,1)=[#V_1;#V_2;#V_3;Nsubdomain]. If using DG, in T we also
%     save the indeces of the edges of the triangles.
% E: edge matrix. The first and second rows contain indices of the starting and ending point of the edge, 
%   the third and fourth rows contain the starting and ending parameter values, 
%   the fifth row contains the edge segment number, and the sixth and seventh row contain 
%   the indeces of subdomains on the left and on the right side.
% Pb: vertex matrix of the degrees of freedom (If using P1, Pb=P).
% Tb: element matrix with indeces degrees of freedom on each element ((If using P1, Tb=T)
% Eb: edge matrix. In Eb the first and second rows contain indices of the starting and
% ending vertices
% the third and fourth rows contains triangle on the left and on the right
% according to anticlockwise orientation


% author: Tommaso Vanzan
%-------------------------------------------------------------------------

function [P,T,Pb,Tb,E,Eb]=generate_mesh_2D(geo,basis_type)
%===================================================
%201: 2D linear nodal FE on triangular mesh
%===== Generation Mesh
if(geo.regular==0)
    [P,E,T]=initmesh(geo.fun,'Hmax',geo.h);%general triangular mesh
else
    [P,E,T] = poimesh(geo.fun,1/geo.h,1/geo.h);%regular triangular mesh
end
    Pb=P;% DOFs correspond to vertex triangles
    Tb=T(1:end-1,:);
    [Eb]=buildmatrixEb_P1_FEM(P,T,E);
end
