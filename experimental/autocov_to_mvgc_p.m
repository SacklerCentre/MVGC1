function F = autocov_to_mvgc_p(G,x,y)

% Calculate partial Granger causality by forward-lagging conditioning variable.

[n,~,p] = size(G);

p = p-1; % we lose the last autocovariance lag

xy = [x y];
z  = 1:n; z(xy) = []; % indices of other variables (to condition out)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward lagging corresponds to a right-shift of the xy,z autocovariance %
% and a left-shift of the z,xy autocovariance.                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GG = G(:,:,1:p);

% note: forward-lagged G_0 is G_{-1} = transpose of G_1'
GG(xy,z,1) = G(z,xy,2)';
GG(z,xy,1) = G(z,xy,2);

% forward-lag remainder of G
for k = 2:p
    GG(xy,z,k) = G(xy,z,k-1);
    GG(z,xy,k) = G(z,xy,k+1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = autocov_to_mvgc(GG,x,y);
