function val=winner(X,Y,M,Xi)
%%WINNER computes the Frobenious inner product weighted by the matrices M and Xi.

val=trace(X'*M*Y*Xi);
end