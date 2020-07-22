
ntrials   = 10;     % number of trials
nobs      = 1000;   % number of observations per trial
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

AT = var9_test;
[nvars,~,p] = size(AT);

% Residuals covariance matrix.

%VT = eye(nvars);
C = tril(randn(nvars));
VT = C*C';
            
X  = var_to_tsdata(AT,VT,nobs,ntrials);                     

[A,V,~,AB,VB] = tsdata_to_var_bootstrap(X,p,nsamps);
