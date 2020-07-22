
ntrials   = 1;      % number of trials
m      = 100000;  % number of observations per trial
nsamps    = 10;     % number of empirical samples

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)

acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)

tstat     = 'chi2'; % statistical test for MVGC: 'FEMP' for Granger's FEMP-test, 'chi2' for Geweke's chi2 test or leave empty for default
alpha     = 0.05;   % significance level for all statistical tests
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

seed      = 1817;   % random seed (0 for unseeded)

x = 5;
y = 3;

%-------------------------------------------------------------------------------

% Seed random number generator.

rng_seed(seed);

% Get VAR coefficients for 5-node test network.

A = var9_test;
[n,~,p] = size(A);  % number of variables

AA = reshape(A,n,n*p);
Y = randn(n,p,m);
YY = reshape(Y,n*p,m);
X = AA*YY + randn(n,m);

%---------------

AA1 = X/YY;  % estimated coeffs
E1 = X - AA1*YY; % residuals

%---------------

EB = E1(:,randi(m,1,m)); % subsampled residuals
XB = X-E1+EB;            % bootstrap
AAB = XB/YY;             % bootstrap coefficients

A1 = reshape(AA1,n,n,p);
AB = reshape(AAB,n,n,p);


