
function C = var2con(A)

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');

C = false(n);
for k = 1:p
    C = C | abs(A(:,:,k)) > eps;
end
