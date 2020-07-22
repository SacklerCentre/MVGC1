function [A,SIG,E] = regress(X,Y,p)

% Regress X on p lags of Y, returning coefficients A, residuals covariance
% matrix SIG and (optionally) residuals E. On return, test for A a NaN
% (bad regression).

[nx,m,N]   = size(X);
[ny,my,Ny] = size(Y);
assert(my == m && Ny == N,'time series don''t match');
assert(p < m,'too many lags');
p1 = p+1;

SIG = NaN;
E   = NaN;

X = demean(X);
Y = demean(Y);

M = N*(m-p);

% stack lags

X0 = reshape(X(:,p1:m,:),nx,M); % concatenate trials for unlagged observations
YL = zeros(ny,p,M);
for k = 1:p
    YL(:,k,:) = reshape(Y(:,p1-k:m-k,:),ny,M); % concatenate trials for k-lagged observations
end
YL = reshape(YL,ny*p,M);        % stack lags

oldstate = warning('off','all'); lastwarn(''); % suppress all warnings
A = X0/YL;                                     % OLS using QR decomposition
[~,warnid] = lastwarn; warning(oldstate);      % check if there were warnings and restore warning state
if ~(isempty(warnid) && all(isfinite(A(:))))   % regerssion failed (something went wrong - we don't really care what... the result is unreliable!)
    A = NaN;
    return;
end

if nargout > 1
    E   = X0-A*YL;              % residuals
    SIG = (E*E')/(M-1);         % residuals covariance matrix (unbiased estimator)
    E   = reshape(E,nx,m-p,N);  % put residuals back into per-trial form
end

A = reshape(A,nx,ny,p);         % so A(:,:,k) is the k-lag coefficients matrix
