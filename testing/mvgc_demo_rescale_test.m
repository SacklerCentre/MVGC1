
ntrials   = 10;     % number of trials
nobs      = 1000;   % number of observations per trial

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'OLS';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 20;     % maximum model order for model order estimation

acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)

tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

fs        = 200;    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)

seed      = 0;      % random seed (0 for unseeded)

%% Generate VAR test data (<mvgc_schema.html#3 |A3|>)
%
% _*Note:*_ This is where you would read in your own time series data; it should
% be assigned to the variable |X| (see below and <mvgchelp.html#4 Common
% variable names and data structures>).

% Seed random number generator.

rng_seed(seed);

% Get VAR coefficients for 5-node test network.

AT = var5_test;
nvars = size(AT,1); % number of variables

% Residuals covariance matrix.

SIGT = eye(nvars);

% Generate multi-trial VAR time series data with normally distributed residuals
% for specified coefficients and covariance matrix.

ptic('\n*** var_to_tsdata... ');
X = var_to_tsdata(AT,SIGT,nobs,ntrials);
ptoc;
X1 = X + mu;
X2 = a*X1;
clear X

%% Model order estimation (<mvgc_schema.html#3 |A2|>)

% Calculate information criteria up to specified maximum model order.

ptic('\n*** tsdata_to_infocrit\n');
[AIC1,BIC1,moAIC1,moBIC1] = tsdata_to_infocrit1(X1,momax,icregmode);
ptoc('*** tsdata_to_infocrit took ');

ptic('\n*** tsdata_to_infocrit\n');
[AIC2,BIC2,moAIC2,moBIC1] = tsdata_to_infocrit1(X2,momax,icregmode);
ptoc('*** tsdata_to_infocrit took ');

maxabs(AIC1-AIC2)
maxabs(BIC1-BIC2)

% Plot information criteria.

figure(1); clf;
plot_tsdata([AIC1 BIC1]',{'AIC','BIC'},1/fs);
title('Model order estimation');
figure(2); clf;
plot_tsdata([AIC2 BIC2]',{'AIC','BIC'},1/fs);
title('Model order estimation');
