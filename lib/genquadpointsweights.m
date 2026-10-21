function data=genquadpointsweights(data,M,Nc,SC)

if SC==1
    [x,ww]=lgwt(Nc,-1,1); % Nodes and weights in 1D using Gauss-Legendre
    w=ww;
    for i=1:M-1
           w=kron(ww,w); %take kronecker product to multiply the weights.
    end
    data.weights=w/(sum(w)); %normalize them
    chi=GenerateTensorNodes(x,M); %take kronecker product of the nodes.
    data.chi=fliplr(chi);
    data.N=size(chi,2);
elseif SC==0
    data.weights=(1/Nc)*ones(Nc,1);
    data.chi=rand(Nc,M);
    data.N=Nc;
end