function X = dlyap_vec(A,Q)

% 'Vectorised' version of dlyap... may not be very efficient.

n = size(A,1);
X = reshape((eye(n*n)-kron(A,A))\Q(:),n,n);
