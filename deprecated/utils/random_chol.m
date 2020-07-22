function C = random_chol(n,rd)

% C = random_chol(n,rmode,rseed)
%
% Return random Cholesky matrix C for n variables
%
% n         Process dimension
% rd        Residuals covariance random distribution

C = tril(ones(n)).*rd.dist(n,n,rd);
