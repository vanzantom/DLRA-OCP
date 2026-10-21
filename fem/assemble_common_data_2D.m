function  [geo]=assemble_common_data_2D(geo,para)
%-------------------------------------------------------------------------
%assemble_mesh_data: assemble mesh-dependent quantities such as the mass
%matrix $Ms$, the discretization of the right hand side $b$ and of the target
%state $y_d$.

% author: Tommaso Vanzan
%-------------------------------------------------------------------------

geo.Nh=length(geo.int);
Ms=sparse(geo.Nv,geo.Nv); % Pre allocat mass matrix
Ms=assemble_matrix_2D(@(x,y,El) 1,geo.P,geo.T,geo.Tb,geo.Tb,geo.Nv,geo.Nv,para.basis_type,0,0,para.basis_type,0,0,para.order_Gauss);%assemble mass matrix 
Ms=Ms(geo.int,geo.int);
geo.Ms=Ms;
end