function points=GenerateTensorNodes(x,M)
%%GENERATETENSORNODES computes the tensor product nodes where x is the set of nodes in 1D, M is number of dimensions
%Input:-array x containing nodes in 1D.
%      -Integer M, being the dimension of the tensor product.
%Output:-array points containing the tensorized nodes

if M<=1
    points=x;
elseif M==2
    [X1,X2]=ndgrid(x);
    points=[X1(:),X2(:)];
elseif M==3
    [X1,X2,X3]=ndgrid(x);
    points=[X1(:),X2(:),X3(:)];
elseif M==4
    [X1,X2,X3,X4]=ndgrid(x);
    points=[X1(:),X2(:),X3(:),X4(:)];
elseif M==5
    [X1,X2,X3,X4,X5]=ndgrid(x);
    points=[X1(:),X2(:),X3(:),X4(:),X5(:)];
elseif M==6
    [X1,X2,X3,X4,X5,X6]=ndgrid(x);
    points=[X1(:),X2(:),X3(:),X4(:),X5(:),X6(:)];
elseif M==7
    [X1,X2,X3,X4,X5,X6,X7]=ndgrid(x);
    points=[X1(:),X2(:),X3(:),X4(:),X5(:),X6(:),X7(:)];
elseif M==8
    [X1,X2,X3,X4,X5,X6,X7,X8]=ndgrid(x);
    points=[X1(:),X2(:),X3(:),X4(:),X5(:),X6(:),X7(:),X8(:)];
elseif M==9
    [X1,X2,X3,X4,X5,X6,X7,X8,X9]=ndgrid(x);
    points=[X1(:),X2(:),X3(:),X4(:),X5(:),X6(:),X7(:),X8(:),X9(:)];
elseif M==10
    [X1,X2,X3,X4,X5,X6,X7,X8,X9,X10]=ndgrid(x);
    points=[X1(:),X2(:),X3(:),X4(:),X5(:),X6(:),X7(:),X8(:),X9(:),X10(:)];
end
end    