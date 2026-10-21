function val=wnorm(X,M,Xi)
%%WNORM computes the weighted Frobenious norm with matrices M and Xi.

val=sqrt(winner(X,X,M,Xi));
end