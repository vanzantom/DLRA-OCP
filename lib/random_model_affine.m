function [data]=random_model_affine(data,geo,para)
%Generates the basis matrices for the random diffusion coefficient

    Asupp=assemble_matrix_2D(@(x,y,El) 1,geo.P,geo.T,geo.Tb,geo.Tb,geo.Nv,geo.Nv,para.basis_type,1,0,para.basis_type,1,0,para.order_Gauss); %Assemble the matrix of the contribution of x gradient
    Asupp=Asupp+assemble_matrix_2D(@(x,y,El) 1,geo.P,geo.T,geo.Tb,geo.Tb,geo.Nv,geo.Nv,para.basis_type,0,1,para.basis_type,0,1,para.order_Gauss);  
    data.Abasis{1}=data.kcoeff*Asupp(geo.int,geo.int);

    l=1;k=1;
    for i=1:data.p       
        Asupp=assemble_matrix_2D(@(x,y,El) (1/(k^2+l^2))*sin(pi*k*x)*sin(pi*l*y),geo.P,geo.T,geo.Tb,geo.Tb,geo.Nv,geo.Nv,para.basis_type,1,0,para.basis_type,1,0,para.order_Gauss); %Assemble the matrix of the contribution of x gradient
        Asupp=Asupp+assemble_matrix_2D(@(x,y,El) (1/(k^2+l^2))*sin(pi*k*x)*sin(pi*l*y),geo.P,geo.T,geo.Tb,geo.Tb,geo.Nv,geo.Nv,para.basis_type,0,1,para.basis_type,0,1,para.order_Gauss);  
        data.Abasis{i+1}=data.kcoeff*Asupp(geo.int,geo.int);
        if(mod(i,2)==1)
            l=l+1;
        else
            k=k+1;
        end
    end
end