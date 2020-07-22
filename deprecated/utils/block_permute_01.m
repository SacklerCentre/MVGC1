%% block_permute
%
% Randomly permute non-overlapping blocks of time series data
%
% <matlab:open('block_permute.m') code>
%
%% Syntax
%
%     XP = block_permute(X,b)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     X          multi-trial time series data
%     b          block size (observations)
%
% _output_
%
%     XP         time series data with randomly permuted blocks
%
%% Description
%
% Randomly permute non-overlapping blocks of size |b| of the input (possible
% multi-trial) time series data |X|. Any short "leftover" block is inserted
% between a randomly chosen pair of blocks.
%
% *_Note:_* For multi-trial time series data, currently each trial is permuted
% independently, although it may make more sense to swap blocks between trials,
% since they are assumed sampled from the same process; this may be implemented
% in a future release.
%
%% See also
%
% <tsdata_to_mvgc_permtest.html |tsdata_to_mvgc_permtest|> |
% <tsdata_to_mvgc_pwc_permtest.html |tsdata_to_mvgc_pwc_permtest|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function XP = block_permute_01(X,b)

[n,m,N] = size(X);

assert(isscalar(b) && isint(b) && b > 0 && b <= m, 'block size must be a positive integer less than the data length');

XP = zeros(n,m,N);
if b == 1
    for r = 1:N
        XP(:,:,r) = X(:,randperm(m),r);   % permute
    end
else
    nb  = floor(m/b);                     % number of blocks
    mm  = nb*b;                           % permutation sequence length
    for r = 1:N
        XX = reshape(X(:,1:mm,r),n,b,nb); % stack first nb blocks...
        XX = XX(:,:,randperm(nb));        % ...and permute them
        XL = X(:,mm+1:m,r);               % short leftover block of length m-mm
        i  = randi([0 nb]);               % random insertion point for leftover block
        % unstack i blocks, insert leftover block and unstack remaining nb-i blocks
        XP(:,:,r) = [reshape(XX(:,:,1:i),n,i*b) XL reshape(XX(:,:,i+1:nb),n,(nb-i)*b)];
    end
end
