function Ynew=realignment(Y)
%%REALIGNMENT receives an array of structures and realignes each component
%%by solving Procrustes problems

Nt=length(Y);
Ynew = repmat(struct('U', [], 'S', [], 'V', []), Nt, 1);
Ynew(1).U=Y(1).U;   Ynew(1).V=Y(1).V;   Ynew(1).S=Y(1).S; 
for i=2:Nt
    AtB=Y(i).U'*Ynew(i-1).U;
    [Ut,~,Vt]=svd(AtB);
    Q1=Ut*Vt';
    Ynew(i).U=Y(i).U*Q1;
    AtB=Y(i).V'*Ynew(i-1).V;
    [Ut,~,Vt]=svd(AtB);
    Q2=Ut*Vt';
    Ynew(i).V=Y(i).V*Q2;
    Ynew(i).S=Q1'*Y(i).S*Q2;
end

