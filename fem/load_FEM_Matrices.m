function [data,pars]=load_FEM_Matrices(filename,data,pars)

load(filename,'Abasis','M','Ml','LM','y0','force');
data.Abasis=Abasis(1:data.p+1);
pars.M=M;
pars.Ml=Ml;
pars.LM=LM;
if size(y0,2)==1        %deterministic initial condition
    data.Y0=repmat(y0,1,data.N);
else
    data.Y0=y0;
end
data.force=force;

end