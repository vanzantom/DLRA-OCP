function [Ar]=GetMatrices(data)
%%GETMATRICES assembles the random matrices
%Input: -  data structure containing the stiffness matrices of the affine
%          expansion
%Output: - Ar structure of length N. Ar{1} is one realization of the
%          stiffness matrix A(\omega).
                            
                
N=data.N;
Ar=cell(N,1);
Abasis=data.Abasis;     %for better parallelization
chi=data.chi;
for i=1:N
       Ar{i}=assemble_affine(Abasis,chi(i,:));
end
