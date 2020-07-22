%% block_permute
%
% Randomly permute non-overlapping blocks of time series data
%
% <matlab:open('block_permute.m') code>
%
%% Syntax
%
%     XP = block_permute(X,bsize)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     X          multi-trial time series data
%     bsize      block size (observations)
%
% _output_
%
%     XP         time series data with randomly permuted blocks
%
%% Description
%
% Randomly permute non-overlapping blocks of size |bsize| of the
% multivariate time series data |X|, which may be multi-trial. The block
% size |bsize| must divide the number of observations in |X| exactly.
%
% For permutation testing of a VAR time series, it is recommended that the
% block size be set to the model order (see e.g. <mvgc_demo_permtest.html
% |mvgc_demo_permtest|>).
%
% *_Note:_* For multi-trial time series data, currently each trial is permuted
% independently, although it may make more sense to swap blocks between trials,
% since they are assumed sampled from the same process; this may be implemented
% in a future release.
%
%% See also
%
% <tsdata_to_var_permtest.html |tsdata_to_var_permtest|> |
% <tsdata_to_mvgc_permtest.html |tsdata_to_mvgc_permtest|> |
% <tsdata_to_mvgc_pwc_permtest.html |tsdata_to_mvgc_pwc_permtest|> |
% <tsdata_to_smvgc_permtest.html |tsdata_to_smvgc_permtest|> |
% <tsdata_to_smvgc_pwc_permtest.html |tsdata_to_smvgc_pwc_permtest|> |
% <mvgc_demo_permtest.html |mvgc_demo_permtest|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function XP = block_permute(X,bsize)

[n,m,N] = size(X);

assert(isscalar(bsize) && isint(bsize) && bsize > 0, 'block size must be a positive integer');

nblocks  = floor(m/bsize); % number of blocks

assert(nblocks*bsize == m,'block size must divide number of observations exactly');

XP = zeros(n,m,N);
for r = 1:N
    XX = reshape(X(:,:,r),n,bsize,nblocks);             % stack blocks
    XP(:,:,r) = reshape(XX(:,:,randperm(nblocks)),n,m); % permute and unstack
end
