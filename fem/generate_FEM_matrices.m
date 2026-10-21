function generate_FEM_matrices(filename,geo,data,paraFEM)

data = random_model_affine(data, geo, paraFEM);   % assemble A basis
M = assemble_matrix_2D(@(x, y, El) 1, geo.P, geo.T, geo.Tb, geo.Tb, geo.Nv, geo.Nv, paraFEM.basis_type, 0, 0, paraFEM.basis_type, 0, 0, paraFEM.order_Gauss);
M = M(geo.int, geo.int);                          % mass matrix
Ml= diag(sum(M, 2));                              % lumped mass matrix
[ LM ] = chol(M);
Abasis=data.Abasis;

                                            % Evaluate initial state
force=zeros(data.Nh,1);


Cov=@(x,y) 0.1*exp(-norm(x-y,2)^2);
y0 = sampleGP2D(geo, Cov, data.N);

for k=1:length(geo.int)
        force(k)=data.f(geo.Pb(1,geo.int(k)),geo.Pb(2,geo.int(k)));
end

save(filename,'Abasis','M','Ml','LM','force','y0');
end