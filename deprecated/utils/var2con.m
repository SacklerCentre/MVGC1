function C = var2con(A)

% C = var2con(A)
%
% Return connectivity matrix (directed graph) for VAR coefficients matrix.
% Self-connections are ignored. Rows index "to" node, columns "from" node.

[n,~,p] = size(A);
C = false(n,n);
for r = 1:p
    C = C | (abs(A(:,:,r)) > eps);
end
C = C+0;            % make numeric
C(1:n+1:n*n) = NaN; % NaNs on diagonal

