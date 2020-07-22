%% autocov_xform
%
% Transform autocovariance sequence for reduced regression
%
%% Syntax
%
%     [GR,AR,SIGR] = autocov_xform(G,x)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     G          autocovariance sequence
%     x          vector of indices of variables for reduced regression
%
% _output_
%
%     GR         transformed autocovariance sequence
%     AR         VAR coefficients matrix for reduced regression
%     SIGR       residuals covariance matrix for reduced regression
%
%% Description
%
% Transforms an autocovariance sequence |G| for new variable defined as residual
% of the reduced regression of the sub-process indexed by |x|. Returns the
% transformed autocovariance sequence |GR| along with the VAR coefficients
% matrix |AR| and residuals covariance matrix |SIGR| of the reduced regression.
%
% Note that the indexing of the returned matrices is _not_ the same as for |G|;
% specifically, indices |1:length(x)| index the variables originally indexed by
% |x| in |G|.
%
% This function is crucial to the calculation of spectral causality in the
% conditional case; see <autocov_to_smvgc.html |autocov_to_smvgc|>,
% <autocov_to_smvgc_pwc.html |autocov_to_smvgc_pwc|>. In theory, if the original
% autocovariance sequence is calculated to |q| lags - under the assumption that
% it may not have decayed sufficiently for |k < q| lags (see
% <var_to_autocov.html |var_to_autocov|>) - then the transformed autocovariance
% sequence should be calculated to |2q| lags. In practice we find that
% calculating to |q| lags is generally sufficient for good accuracy. To
% calculate |GR| to higher lags, the simplest option is to reduce the |acdectol|
% parameter in the call to <var_to_autocov.html |var_to_autocov|> (e.g. squaring
% it will effectively double the number of lags |q| to which |G| and hence |GR|
% is calculated).
%
%% References
%
% [1] L. Barnett and A. K. Seth, <matlab:openpdf('mvgc_doc.pdf') The MVGC
% Multivariate Granger Causality Toolbox>, _in preparation_, Aug. 2012.
%
%% See also
%
% <var_to_autocov.html |var_to_autocov|> |
% <autocov_to_var.html |autocov_to_var|> |
% <autocov_to_smvgc.html |autocov_to_smvgc|> |
% <autocov_to_smvgc_pwc.html |autocov_to_smvgc_pwc|>
%
%% Copyright notice
%
% _(C) Lionel Barnett, 2012. See file <matlab:edit('license.txt')
% license.txt> in root directory for licensing terms._
%
%%

function [GR,AR,SIGR] = autocov_xform_old(G,x)

[n,~,q1] = size(G);

x = x(:)';          % vectorise
y = 1:n; y(x) = []; % indices of other variables

assert(all(x >=1 & x <= n),'some indices out of range');
assert(~isempty(y),'no variables left!');

q = q1-1;

% transform autocov by reduced regression

xy = [x y];
G = G(xy,xy,:); % rearrange
nx = length(x);
ny = length(y);
x = 1:nx;
y = nx+1:n;

[AR,SIGR] = autocov_to_var(G(x,x,:)); % reduced regression
AR = reshape(cat(3,eye(nx),-AR),nx,q1*nx);

GF = reshape(G(y,x,:),ny,q1*nx);                             % forward  autocovariance sequence
GB = reshape(permute(flipdim(G(x,y,:),3),[1 3 2]),q1*nx,ny); % backward autocovariance sequence

GR = zeros(n,n,q1);   % transformed autocovariance sequence
GR(x,x,1) = SIGR;     % just the reduced residuals covariance, since residuals serially uncorrelated
GR(y,y,:) = G(y,y,:); % just the y autocovariance
for k = 0:q-1
    GR(x,y,k+1) = AR(:,1:(k+1)*nx)*GB((q-k)*nx+1:q1*nx,:) + AR(:,(k+1)*nx+1:q1*nx)*GF(:,nx+1:(q1-k)*nx)';
end
GR(x,y,q1) = AR*GB;
for k = 0:q
    GR(y,x,k+1) = GF(:,k*nx+1:q1*nx)*AR(:,1:(q1-k)*nx)';
end
