function X = v_var_to_tsdata_old(A,SIG,N,ms)

% X = v_var_to_tsdata(A,SIG,N,ms)
%
% Return N independent random Gaussian VAR time series of length m with
% coefficients A and covariance matrix SIG.
%
% First form: if ms >= 1 it is taken to be the the number of steps 'ms'
% to truncate until (assumed) stationarity.
%
% Second form: if ms < 1 it is assumed to be the spectral radius 'rho'
% of A and the number of steps to stationarity is (over)estimated
% automatically; set 'decayfac' larger for longer settle time.

% test:
%
% n=2; p=1; m = 1000; A = zeros(n,n,p,m);
% for t = 1:m, A(:,:,:,t) = [0.3 sin(2*pi*(t/1000-1)); 0.2 0]; end

if nargin < 4 || isempty(ms), ms = 0; end % default: no truncation

[n,n1,p,m] = size(A);
assert(n1 == n);
[n2,n3,m1] = size(SIG);
assert(n2 == n && n3 == n);
statsig = m1 == 1;
if statsig
    [SIG,cholp] = chol(SIG,'lower');
    assert(cholp == 0,'covariance matrix not positive-definite');
else
    assert(m1 == m);
end

X = zeros(n,m-ms,N);

for r = 1:N % for each realization

    X(:,:,r) = v_genvar(A,SIG,n,p,m,ms,statsig);

end

function X = v_genvar(A,SIG,n,p,m,ms,statsig)

% initialise to Gaussian white noise

if statsig
    X = SIG*randn(n,m); % "SIG" is actually Cholesky matrix
else
    X = zeros(n,m);
    for t = 1:m % for each time step
        [C,cholp] = chol(SIG(:,:,t),'lower');
        assert(cholp == 0,'covariance matrix not positive-definite at time step %d',t);
        X(:,t) = C*randn(n,1);
    end
end

% loop through time steps

for t = p+1:m % for each time step

    for k = 1:p % for each lag
        X(:,t) = X(:,t) + A(:,:,k,t)*X(:,t-k);
    end

end

% truncate first ms

X = X(:,ms+1:m);
