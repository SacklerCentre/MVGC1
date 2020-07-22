%% autocov_xform
%
% Transform autocovariance sequence for reduced regression
%
%% Syntax
%
%     GR = autocov_xform(G,AR,SIGR)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     G          autocovariance sequence
%     AR         VAR coefficients matrix for reduced regression
%     SIGR       residuals covariance matrix for reduced regression
%
% _output_
%
%     GR         transformed autocovariance sequence
%
%% Description
%
% Returns the autocovariance sequence |GR| for a new variable defined as the
% residuals of a reduced regression, for a VAR with autocovariance sequence |G|.
% |AR| and |SIGR| are the coefficients matrices and residuals covariance matrix
% respectively of the reduced regression; the reduced regression is assumed to
% correspond to the first |size(AR,1)| indices of |G|.
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
% <autocov_to_smvgc_pwc.html |autocov_to_smvgc_pwc|> |
% <cpsd_xform.html |cpsd_xform|>
%
%% Copyright notice
%
% _(C) Lionel Barnett, 2012. See file <matlab:edit('license.txt')
% license.txt> in root directory for licensing terms._
%
%%

function GR = autocov_xform_old2(G,AR,SIGR)

[n,~,q1] = size(G);
q = q1-1;

[nx,nx1,~] = size(AR);
assert(nx1 == nx,'reduced VAR coefficients matrix has bad shape');
assert(nx <= n,'reduced VAR coefficients matrix appears to be for more variables than autocovariance sequence');

[n1,n2] = size(SIGR);
assert(n1 == n2,'reduced VAR residuals covariance matrix not square');
assert(n1 == nx ,'reduced VAR residuals covariance matrix doesn''t match reduced VAR coefficients matrix');

ny = n-nx;
x = 1:nx;
y = nx+1:n;

% transform autocov by reduced regression

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
