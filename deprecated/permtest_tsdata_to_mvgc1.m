function FP = permtest_tsdata_to_mvgc1(U,x,y,p,bsize,nsamps,regmode,acmaxlags,acdectol)

if nargin < 7, regmode   = []; end % ensure default
if nargin < 8, acmaxlags = []; end % ensure default
if nargin < 9, acdectol  = []; end % ensure default

[n,m,N] = size(U);
assert(m > p,'too many lags');
p1 = p+1;
M = N*(m-p);
np = n*p;

x = x(:)'; % vectorise
y = y(:)'; % vectorise
assert(all(x >=1 & x <= n),     'some x indices out of range');
assert(all(y >=1 & y <= n),     'some y indices out of range');
assert(isempty(intersect(x,y)), 'x and y indices must be distinct');

w = 1:n; w(x) = [];     % w = yz
nx = length(x);
ny = length(y);
nw = length(w);

assert(isscalar(bsize) && isint(bsize) && bsize > 0, 'block size must be a positive integer');
nblocks  = floor(m/bsize); % number of blocks
assert(nblocks*bsize == m,'block size must divide number of observations exactly');

FP = nan(nsamps,1);

EP = zeros(n,M);
AA = zeros(n,np);

UL = zeros(n,p,M);
for k = 1:p
    UL(:,k,:) = reshape(U(:,p1-k:m-k,:),n,M); % concatenate trials for k-lagged observations
end
UL = reshape(UL,np,M);          % stack lags
W0 = reshape(U(w,p1:m,:),nw,M); % concatenate trials for unlagged observations
AA(w,:) = W0/UL;                % OLS using QR decomposition
EP(w,:) = W0-AA(w,:)*UL;        % residuals

UP = U; % with permuted y-variable

for s = 1:nsamps
    fprintf('MVGC: permutation test sample %d of %d',s,nsamps);
    
    % generate permutation time series: "block permute" the y-variable per-trial
    
    for r = 1:N
        Yr = reshape(U(y,:,r),ny,bsize,nblocks);            % stack blocks
        UP(y,:,r) = reshape(Yr(:,:,randperm(nblocks)),ny,m); % permute blocks and unstack
    end
    
    UPL = zeros(n,p,M);
    for k = 1:p
        UPL(:,k,:) = reshape(UP(:,p1-k:m-k,:),n,M); % concatenate trials for k-lagged observations
    end
    UPL = reshape(UPL,np,M);        % stack lags
    X0 = reshape(U(x,p1:m,:),nx,M); % concatenate trials for unlagged observations
    AA(x,:) = X0/UPL;               % OLS using QR decomposition
    EP(x,:) = X0-AA(x,:)*UPL;       % residuals
    
    if isbad(AA), fprintf(' *** VAR estimation failed\n'); continue; end % something went badly wrong

    AP = reshape(AA,n,n,p);         % so AP(:,:,k) is the k-lag coefficients matrix
    SIGP = (EP*EP')/(M-1);          % residuals covariance matrix
    
    % calculate permutation test MVGC
    
    [G,res] = var_to_autocov(AP,SIGP,acmaxlags,acdectol);
    if res.error, fprintf(' *** bad VAR: %s\n',res.errmsg); continue; end
    FP(s) = autocov_to_mvgc(G,x,y);
    
    fprintf('\n');
end
