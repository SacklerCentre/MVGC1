function [AB,SIGB] = tsdata_to_var_bootstrap_01(X,p,nsamps)

[n,m,N] = size(X);
assert(p < m,'too many lags');
p1 = p+1;
M = N*(m-p);
np = n*p;

AB   = nan(n,n,p,nsamps);
SIGB = nan(n,n,nsamps);

X = demean(X); % no constant term

% estimate VAR coefficients

X0 = reshape(X(:,p1:m,:),n,M); % concatenate trials for unlagged observations
XL = zeros(n,p,M);
for k = 1:p
    XL(:,k,:) = reshape(X(:,p1-k:m-k,:),n,M); % concatenate trials for k-lagged observations
end
XL = reshape(XL,np,M);         % stack lags
A = X0/XL;                     % OLS using QR decomposition
if isbad(A), return; end       % something went badly wrong

% calculate predictions and residuals

m   = m-p;              % we lose p observations
XP  = A*XL;             % predictions
E   = X0-XP;            % residuals: so X0 = XP + E
E   = reshape(E,n,m,N); % put residuals back into per-trial form

% bootstrap VAR parameters

for s = 1:nsamps
    fprintf('estimating bootstrap parameters %d of %d',s,nsamps);
    XB = XP + reshape(E(:,randi(m,1,m),:),n,M); % bootstrap: subsample residuals with replacement and add to predictions
    A1 = XB/XL;                                 % OLS using QR decomposition
    if isbad(A1), fprintf(' *** VAR estimation failed\n'); continue; end   % something went badly wrong
    AB(:,:,:,s) = reshape(A1,n,n,p);            % bootstrap VAR coefficients
    E1 = XB-A1*XL;                              % bootstrap residuals
    SIGB(:,:,s) = (E1*E1')/(M-1);               % bootstrap residuals covariance matrix
    fprintf('\n');
end
