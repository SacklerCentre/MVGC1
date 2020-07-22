function [A,SIG,E] = lregress(Y,X,p)

% regress Y against p lags of X

[nx,m]  = size(X);
[ny,m1] = size(Y);
assert(m1 == m);

M = m-p;

% stack lags

Y = Y(:,p+1:m);
XL = zeros(nx,p,M);
for k = 1:p
    XL(:,k,:) = X(:,p+1-k:m-k);
end
XL = reshape(XL,nx*p,M);

A = Y/XL;

if nargout > 1
    E   = Y-A*XL;
    SIG = (E*E')/(M-1);
end

A = reshape(A,ny,nx,p);
