function A = lesion(A,C)

% A = lesion(A,C)
%
% Lesion connections in list C (format: rows of [to,from]) in VAR
% coefficients matrix A.

if isempty(C)
    return;
end

p = size(A,3);
assert(size(C,2) == 2,'connections matrix must have two columns');

nc = size(C,1);
for c = 1:nc
    for k = 1:p
        A(C(c,1),C(c,2),k) = 0;
    end
end

for k = p:-1:2
    AA = A(:,:,k);
    if any(AA(:)), break; end
end
AA = A(:,:,k);
if k == 2 && ~any(AA(:))
    fprintf(2,'WARNING: no connections left!\n');
    A = [];
elseif k < p
    A = A(:,:,1:k);
end

