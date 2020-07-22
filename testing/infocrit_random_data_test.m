nvars    = 8;    % number of variables
nobs     = 1000; % number of observations
ntrials  = 1;    % number of trials

maxorder = 30;   % maximum model order to test

rseed    = 0;     % random seed - zero for unseeded)

regmode  = 'LWR'; % regression mode ('LWR' or 'OLS')

%-------------------------------------------------------------------------------

% Generate random Gaussian time series (white noise)

oseed = rng_seed(rseed);
X = randn(nvars,nobs,ntrials);
rng_restore(oseed);

% Calculate AIC, BIC

[aic,bic,moaic,mobic] = tsdata_to_infocrit(X,maxorder,regmode);

% Display

figure(1); clf;
plot_tsdata([aic bic]',{sprintf('AIC (opt = %d)',moaic),sprintf('BIC (opt = %d)',mobic)});
title(sprintf('Model order estimation (mode: %s)',regmode));
