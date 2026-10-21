function val=Bochnorm(dz,Ms,data)
%BOCHNORM computes the square Bochner norm over the time interval (0,T) induced by the Ms-norm.
%It receives the vector dz, the matrix Ms and the parameters Nt and dt in the structure data.

Nt=data.Nt;dt=data.dt;

val=0;
for i=1:Nt
    val=val+dt*(dz(:,i))'*Ms*dz(:,i);
end
