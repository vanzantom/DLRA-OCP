function G = sampleGP2D(geo, C, N)
% SAMPLEGP2D Generate a set of N Gaussian processes realizations on a 2D grid
% and enforce homogeneous Dirichlet boundary conditions.
%
% Inputs:
%   geo : mesh data structuer
%   C : covariance function handle  C(x1,x2)
%   N : number of samples

% Output:
%   G : n-by-n matrix containing the sampled field
   x=geo.P(1,:);
   y=geo.P(2,:);
   n=length(x);
   K = zeros(n,n);
    for i=1:n                   %exploit symmetry
        for j=i:n
            K(i,j)=C([x(i);y(i)],[x(j);y(j)]);
            K(j,i) = K(i,j);
        end
    end
    K=K+1e-10*eye(size(K));
    L = chol(K,'lower');
    G=zeros(geo.Nh,N);
    for i=1:N
        g = L*randn(n,1);           %sample GP
        g = x' .* (1-x)' .* y' .* (1-y)' .* g; %enforce boundary
        G(:,i)=g(geo.int);
    end

end