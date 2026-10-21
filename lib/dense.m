function Y=dense(Ylr)
%%DENSE receives a structure with fields U (mxr),S (rxr),V(nxr) representing a low-rank
%matrix and returns the full matrix m\times n

Y=Ylr.U*Ylr.S*Ylr.V';

end