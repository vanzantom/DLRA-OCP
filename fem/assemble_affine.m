function A=assemble_affine(Abasis,chi)

A=Abasis{1};
for ell=1:length(chi)
    A=A+chi(ell)*Abasis{ell+1};
end