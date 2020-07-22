function GD = autocov_self(G)

% GD = autocov_self(G)
%
% Return matrix of self-autocovariances (without cross-autocovariances).

[n,~,q1] = size(G);

GD = zeros(n,q1);
for k = 1:q1
    GD(:,k) = diag(G(:,:,k));
end
