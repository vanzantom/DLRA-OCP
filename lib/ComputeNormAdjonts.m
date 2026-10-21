function val=ComputeNormAdjonts(P,Q,Ms,Xi,data)
%COMPUTENORMADJOINTS computes the square Bochner norm over the time interval (0,T) induced by the (M,\Xi)-weigthed Frobenious inner product.
%It receives the matrix P (expected to be the exact adjiont) and
%the structure Q (either of fields U,S,V or Ubase,Vbase,dS,dU,dV),the matrix Ms and the parameters Nt and dt in the structure data.

Nt=data.Nt;dt=data.dt;

val=0;
if numel(fieldnames(Q))==3
    for i=1:Nt
        Qdense=dense(Q(i));
        diff=P(:,:,i)-Qdense;
        val=val+dt*winner(diff,diff,Ms,Xi);
    end
elseif numel(fieldnames(Q))==5
    for i=1:Nt
        Qdense=tantodens(Q(i));
        diff=P(:,:,i)-Qdense;
        val=val+dt*winner(diff,diff,Ms,Xi);
    end
else
    fprintf("Adjoint variable is not in the correct format\n");
end
