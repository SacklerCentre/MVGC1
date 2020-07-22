function A = random_var(G,p,wd,wds)

% A = random_coeffs(G,p,wd,wds)
%
% Generate VAR coefficients on an adjacency matrix from a specified random
% distribution. Optionally add in random self-activations.
%
% G          graph adjacency matrix
% p          model order
% wd         weight distribution structure (non-self)
% wds        weight distribution structure (self)
%
% Weight distribution structure:
%
% wd.dist    randomisation function (like randn, etc)
% wd.dfac    decay factor
% wd.???     other parameters (which will depend on particular wd.dist)

n = size(G,1);
assert(size(G,2) == n, 'adjacency matrix not square');

k = nnz(G); % number of edges
W = wd.dist(p,k,wd);
A = zeros(n,p*n);
d = 1;
for r=1:p
    V = zeros(1,n*n);
    V(G) = d*W(r,:);
    A(:,(r-1)*n+1:r*n) = reshape(V,n,n);
    d = wd.dfac*d;
end

if nargin > 3 && ~isempty(wds)
    W = wds.dist(p,n,wds);
    d = 1;
    dgp = 1:n+1:n*n; % diagonal
    for r=1:p
        A(dgp) = d*W(r,:);
        dgp = dgp+n*n;
        d = wds.dfac*d;
    end
end

A = reshape(A,n,n,p);
