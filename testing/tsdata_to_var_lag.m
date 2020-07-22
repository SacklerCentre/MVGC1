function [A,SIG,E] = tsdata_to_var_lag(X,p,nlags)

[n,m,N] = size(X);

P = p*nlags;

assert(P < m,'too many lags');
p1 = p+1;

A   = NaN; % assure a "bad" return value if anything goes wrong (see routine 'isbad')
SIG = NaN; % assure a "bad" return value if anything goes wrong (see routine 'isbad')
E   = NaN; % assure a "bad" return value if anything goes wrong (see routine 'isbad')

X = demean(X); % no constant term

M = N*(m-p);
np = n*p;

% stack lags

X0 = reshape(X(:,p1:m,:),n,M); % concatenate trials for unlagged observations
XL = zeros(n,p,M);
for k = 1:p
    XL(:,k,:) = reshape(X(:,p1-k:nlags:m-k,:),n,M); % concatenate trials for k-lagged observations
end
XL = reshape(XL,np,M);         % stack lags

A = X0/XL;                     % OLS using QR decomposition
if isbad(A); return; end       % something went badly wrong

if nargout > 1
    E   = X0-A*XL;             % residuals
    SIG = (E*E')/(M-1);        % residuals covariance matrix
    E   = reshape(E,n,m-p,N);  % put residuals back into per-trial form
end

A = reshape(A,n,n,p);          % so A(:,:,k) is the k-lag coefficients matrix
