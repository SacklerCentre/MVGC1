function [R,p] = rchol(A)

% "Reverse Cholesky": R = rchol(A) produces an upper triangular matrix R from
% the diagonal and upper triangle of matrix A, satisfying the equation R*R'= A.

[L,p] = chol(rot90(A,2));
if p == 0
    R = rot90(L,2)';
end
