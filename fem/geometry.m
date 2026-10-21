function [geo]=geometry(h,data,para)

geo=create_geometry_2D(h,para.basis_type);  
geo.boundary=generate_boundaryedges(geo.Tb,geo.Eb,geo.E,data.label,para.basis_type);
geo=treat_boundary(geo);
geo=assemble_common_data_2D(geo,para);

end