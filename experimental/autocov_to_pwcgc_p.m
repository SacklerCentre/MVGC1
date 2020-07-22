function F = autocov_to_pwcgc_p(G)

% Calculate pairwise-conditional partial Granger causalities by forward-lagging
% conditioning variables.
%
% Note that we have to forward-lag the autocovariance _separately_ for each pair
% of variables (see mvgc_to_autocov_p). This means that effectively we call the
% regression steps (autocov_to_var) 2*n*(n-1) times. This routine is thus _much_
% less efficient than autocov_to_pwcgc, which only requires n calls. This is
% unavoidable.

[n,~,p] = size(G);

p = p-1; % we lose the last autocovariance lag

F = nan(n);

for j = 1:n
    for i=1:n;
        if i == j, continue; end
        
        ij  = [i j];
        oij = 1:n; oij(ij) = []; % omit i, j (to condition out)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Forward lagging corresponds to a right-shift of the ij,oij autocovariance %
        % and a left-shift of the oij,ij autocovariance.                            % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        GG = G(:,:,1:p);

        % note: forward-lagged G_0 is G_{-1} = transpose of G_1'
        GG(ij,oij,1) = G(oij,ij,2)';
        GG(oij,ij,1) = G(oij,ij,2);

        % forward-lag remainder of G
        for k = 2:p
            GG(ij,oij,k) = G(ij,oij,k-1);
            GG(oij,ij,k) = G(oij,ij,k+1);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        F(i,j) = autocov_to_mvgc(GG,i,j);
    end
end
