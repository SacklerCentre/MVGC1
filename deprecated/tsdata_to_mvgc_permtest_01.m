%% tsdata_to_mvgc_permtest
%
% Return permutation samples of conditional time-domain MVGC (multivariate Granger causality)
%
% <matlab:open('tsdata_to_mvgc_permtest.m') code>
%
%% Syntax
%
%     F = tsdata_to_mvgc_permtest(U,x,y,p,nperms,bsize,regmode,acmaxlags,acdectol)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     U          multi-trial time series data
%     x          vector of indices of target (causee) multi-variable
%     y          vector of indices of source (causal) multi-variable
%     p          model order (number of lags)
%     nperms     number of permutation samples
%     bsize      block size (default: model order)
%     regmode    regression mode (default: as for 'tsdata_to_var')
%     acmaxlags  maximum autocovariance lags to calculate (default: as for 'var_to_autocov')
%     acdectol   autocovariance decay (default: as for 'var_to_autocov')
%
% _output_
%
%     F          vector of permutation Granger causalities
%
%% Description
%
% Return |nperms| permutation samples of MVGC of variable specified by indices
% |y| to variable specified by indices |x|, conditional on all other variables
% in the time series data in |U|. |p| specifies the VAR model order.
%
% Prior to each sample MVGC calculation, non-overlapping blocks of size |bsize|
% of the predictor variable are permuted (see <block_permute.html
% |block_permute|>). The idea is to disrupt causal correlation while retaining
% non-causal correlation structure in the data. The default is to set the block
% size to the model order |p|, since this might be expected to best preserve
% non-causal correlation structure. The regression mode parameter |regmode| is
% passed through to <tsdata_to_var.html tsdata_to_var>, while the parameters
% |acmaxlags| and |acdectol| are passed through to <var_to_autocov.html
% var_to_autocov>; see <autocov_to_mvgc.html autocov_to_mvgc> for MVGC
% calculation.
%
% This routine is designed to generate an empirical null MVGC distribution that
% may be used for permutation testing; see <empirical_pval.html empirical_pval>
% and <mvgc_demo_permtest.html mvgc_demo_permtest>.
%
%% References
%
% [1] L. Barnett and A. K. Seth, <matlab:open('mvgc_doc.pdf') The MVGC
% Multivariate Granger Causality Toolbox>, _in preparation_, Aug. 2012.
%
%% See also
%
% <block_permute.html block_permute> |
% <tsdata_to_var.html tsdata_to_var> |
% <var_to_autocov.html var_to_autocov> |
% <autocov_to_mvgc.html autocov_to_mvgc> |
% <empirical_pval.html empirical_pval> |
% <mvgc_demo_permtest.html mvgc_demo_permtest>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function F = tsdata_to_mvgc_permtest_01(U,x,y,p,nperms,bsize,regmode,acmaxlags,acdectol)

if nargin < 6 || isempty(bsize), bsize = p; end % default is to use model order

if nargin < 7, regmode   = []; end % ensure default
if nargin < 8, acmaxlags = []; end % ensure default
if nargin < 9, acdectol  = []; end % ensure default

[n,m,N] = size(U);
assert(m > p,'too many lags');

x = x(:)'; % vectorise
y = y(:)'; % vectorise
assert(all(x >=1 & x <= n),     'some x indices out of range');
assert(all(y >=1 & y <= n),     'some y indices out of range');
assert(isempty(intersect(x,y)), 'x and y indices must be distinct');

UP = U;
F = nan(nperms,1);
for t = 1:nperms
    fprintf('permutation %d of %d',t,nperms);
    UP(y,:,:) = block_permute(U(y,:,:),bsize);
    [A,SIG] = tsdata_to_var(UP,p,regmode);
    if isbad(A), fprintf(' *** VAR estimation failed\n'); continue; end
    [G,res] = var_to_autocov(A,SIG,acmaxlags,acdectol);
    if res.error, fprintf(' *** bad VAR: %s\n',reg.errmsg); continue; end
    F(t) = autocov_to_mvgc(G,x,y);
    fprintf('\n');
end
