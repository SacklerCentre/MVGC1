function A = lesionall(A,to,from)

% A = lesionall(A,to,from)
%
% Lesion all connections connections to nodes in 'to' list from nodes in
% 'from' list, for VAR coefficients matrix A.

p = size(A,3);
assert(isvector(to) && isvector(from),'node lists must be vectors');

for f = from
    for t = to
        for k = 1:p
            A(t,f,k) = 0;
        end
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

