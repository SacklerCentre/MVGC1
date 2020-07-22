%% MVGC for tapping study

%% Parameters

% ntrials = 1;      % number of trials                 - LB: WE'LL READ THIS IN LATER (LINE 50)
% nobs    = 124;    % number of observations per trial - LB: WE'LL READ THIS IN LATER (LINE 50)
regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)
morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 20;     % maximum model order for model order estimation
acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)
tstat     = 'F';    % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')
fs        = 1;      % sample rate (Hz) - LB: SHOULD BE THE (AVERAGE?) TAPPING
fres      = [];     % frequency resolution (empty for automatic calculation)
% seed    = 0;      % random seed (0 for unseeded)       - LB: DON'T NEED THIS

% LB: MORE FLEXIBLE TO READ DATA FROM A FILE - SEE MATLAB 'readmatrix' FUNCTION
%
% NOTE: DATA GOES IN 'X' VARIABLE, 1 ROW PER TIME SERIES

X(1,:) = [-0.0400000000000000,0.00199999999999978,-0.00800000000000001,-0.0179999999999998,0.00200000000000067,0.0559999999999992,-0.0289999999999999,0.00300000000000011,0.00200000000000067,0.00300000000000011,0.0129999999999999,0.0350000000000001,-0.0399999999999992,-0.0190000000000019,-0.00799999999999912,-0.00799999999999912,-0.0180000000000007,0.0449999999999982,-0.0289999999999999,0.0130000000000017,-0.0190000000000019,-0.00799999999999912,0.0139999999999993,0.00200000000000244,0.00300000000000011,-0.00800000000000267,0.00200000000000244,0.00300000000000011,0.00299999999999656,-0.0399999999999992,0.00200000000000244,0.0240000000000009,-0.0180000000000007,0.0129999999999981,-0.0399999999999992,0.00300000000000011,0.0129999999999981,-0.0180000000000007,-0.0730000000000004,0.0140000000000029,0.0450000000000017,0.0129999999999981,-0.0499999999999972,0.0449999999999946,0.00300000000000011,0.00300000000000011,-0.00900000000000034,-0.0289999999999964,0.0129999999999981,-0.0510000000000019,0.0250000000000057,0.00199999999999534,-0.00799999999999557,0.0139999999999958,0.0129999999999981,-0.0189999999999984,-0.0180000000000007,-0.0300000000000011,-0.0399999999999992,0.0450000000000017,-0.00800000000000267,-0.0180000000000007,0.0130000000000052,-0.00799999999999557,0.0240000000000009,-0.0290000000000106,0.0240000000000009,-0.0720000000000027,-0.0189999999999912,0.00300000000000011,-0.0400000000000063,-0.0619999999999976,0.00300000000000011,-0.00799999999999557,-0.0190000000000055,-0.0399999999999920,-0.00700000000000500,-0.0409999999999968,-0.0290000000000106,-0.0189999999999912,0.0129999999999910,0.00300000000000011,-0.00799999999999557,-0.0609999999999928,-0.0290000000000106,0.00200000000000955,-0.00800000000000978,-0.0289999999999964,-0.0510000000000019,-0.00799999999999557,-0.0400000000000063,-0.0289999999999964,-0.0510000000000019,0.00300000000000011,-0.0510000000000019,-0.125000000000000,-98.9750000000000,97.9450000000000,0.0350000000000108,-0.0190000000000055,0.0250000000000057,0.00199999999999534,-0.0390000000000015,0.00199999999999534,-0.0399999999999920,-0.00799999999999557,0.0339999999999918,-0.0399999999999920,-0.0610000000000071,0.0880000000000081,0.0669999999999931,-0.0300000000000011,-0.00799999999999557,-0.00800000000000978,-0.0499999999999972,0.00200000000000955,0.0449999999999875,0.0350000000000108,-0.0180000000000007,-0.0510000000000019,-0.00799999999999557,0.0449999999999875,0.0560000000000116,0.0129999999999910];
X(2,:) = [0.0110000000000001,-0.0620000000000003,0.0240000000000000,-0.0609999999999999,0.0449999999999999,0.00300000000000100,0.0129999999999999,-0.0190000000000001,0.00300000000000011,0.0129999999999999,-0.0180000000000007,-0.00799999999999912,0.0129999999999999,-0.0190000000000019,0.0240000000000009,0.00400000000000134,-0.0300000000000011,0.00199999999999889,-0.0289999999999999,0.0240000000000009,-0.0289999999999999,0.0130000000000017,0.00199999999999889,0.00300000000000011,0.0240000000000009,-0.0190000000000019,0.0240000000000009,-0.0289999999999999,0.00300000000000011,-0.00799999999999912,0.0129999999999981,0.0129999999999981,-0.0289999999999964,-0.00800000000000267,-0.0289999999999964,-0.0190000000000055,0.0130000000000052,-0.0180000000000007,-0.00800000000000267,-0.0189999999999984,0.0240000000000009,-0.0399999999999992,0.0240000000000009,0.0449999999999946,0.00300000000000011,0.0350000000000037,-0.0300000000000011,0.0240000000000009,0.0240000000000009,-0.0190000000000055,-0.0179999999999936,-0.0400000000000063,-0.00799999999999557,0.00300000000000011,-0.00800000000000267,0.0129999999999981,-0.0289999999999964,-0.0300000000000011,-0.00800000000000267,0.00300000000000011,-0.0289999999999964,0.00199999999999534,0.0240000000000009,0.0130000000000052,-0.00799999999999557,-0.0400000000000063,0.00300000000000011,-0.0510000000000019,0.0240000000000009,-0.0400000000000063,-0.0289999999999964,-0.0189999999999912,0.0559999999999974,-0.0820000000000078,-0.0510000000000019,0.0240000000000009,-0.0289999999999964,-0.0399999999999920,0.0129999999999910,-0.0289999999999964,-0.0190000000000055,-0.0609999999999928,-0.00799999999999557,-0.00900000000000034,-0.0180000000000007,0.0129999999999910,-0.0289999999999964,-0.0610000000000071,-0.0189999999999912,0.00199999999999534,0.0240000000000009,-0.0609999999999928,-0.103999999999999,-0.00800000000000978,-0.0399999999999920,-0.0930000000000035,-0.0190000000000055,-0.0399999999999920,-0.0190000000000055,0.0250000000000057,-0.0300000000000011,0.0240000000000009,0.00300000000000011,-0.0300000000000011,0.00300000000000011,-0.00799999999999557,-0.0290000000000106,-0.0300000000000011,0.0140000000000100,-0.00800000000000978,0.0560000000000116,0.0339999999999918,-0.0499999999999972,-0.0300000000000011,0.0139999999999958,0.00200000000000955,-0.00800000000000978,0.00300000000000011,0.0130000000000052,0.0349999999999966,-0.0189999999999912,0.0349999999999966,-0.0400000000000063,0.0559999999999974];

% LB: DON'T NEED ANY OF THIS - IT'S FOR TEST DATA, YOU'VE GOT THE REAL THING :-)
%
% Seed random number generator.
%
% [seed]
%
% Get VAR coefficients for 5-node test network.
%
% AT = var5_test;
% nvars = size(AT,1); % number of variables
%
% Residuals covariance matrix.
%
% SIGT = eye(nvars);
%
% Generate multi-trial VAR time series data with normally distributed residuals
% for specified coefficients and covariance matrix.
%
% ptic('\n*** var_to_tsdata... ');
% X = var_to_tsdata(AT,SIGT,nobs,ntrials);
% ptoc;

% LB: DATA DIMENSIONS SHOULD BE SPECIFIED FROM ACTUAL DATA - LIKE THIS

[nvars,nobs,ntrials] = size(X);

%% Model order estimation (<mvgc_schema.html#3 |A2|>)

% Calculate information criteria up to specified maximum model order.

ptic('\n*** tsdata_to_infocrit\n');
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
ptoc('*** tsdata_to_infocrit took ');

% Plot information criteria.

figure(1); clf;
plot_tsdata([AIC BIC]',{'AIC','BIC'},1/fs);
title('Model order estimation');

% amo = size(AT,3); % actual model order - LB: NOT RELEVANT... WE DON'T KNOW THE "ACTUAL" MODEL ORDER!

fprintf('\nbest model order (AIC) = %d\n',moAIC);
fprintf('best model order (BIC) = %d\n',moBIC);
% fprintf('actual model order     = %d\n',amo);

% Select model order.

if strcmpi(morder,'AIC')
    morder = moAIC;
    fprintf('\nusing AIC best model order = %d\n',morder);
elseif strcmpi(morder,'BIC')
    morder = moBIC;
    fprintf('\nusing BIC best model order = %d\n',morder);
else
    fprintf('\nusing specified model order = %d\n',morder);
end

%% VAR model estimation (<mvgc_schema.html#3 |A2|>)

% Estimate VAR model of selected order from data.

ptic('\n*** tsdata_to_var... ');
[A,SIG] = tsdata_to_var(X,morder,regmode);
ptoc;

% Check for failed regression

assert(~isbad(A),'VAR estimation failed');

% NOTE: at this point we have a model and are finished with the data! - all
% subsequent calculations work from the estimated VAR parameters A and SIG.

%% Autocovariance calculation (<mvgc_schema.html#3 |A5|>)

% The autocovariance sequence drives many Granger causality calculations (see
% next section). Now we calculate the autocovariance sequence G according to the
% VAR model, to as many lags as it takes to decay to below the numerical
% tolerance level, or to acmaxlags lags if specified (i.e. non-empty).

ptic('*** var_to_autocov... ');
[G,info] = var_to_autocov(A,SIG,acmaxlags);
ptoc;

% The above routine does a LOT of error checking and issues useful diagnostics.
% If there are problems with your data (e.g. non-stationarity, colinearity,
% etc.) there's a good chance it'll show up at this point - and the diagnostics
% may supply useful information as to what went wrong. It is thus essential to
% report and check for errors here.

var_info(info,true); % report results (and bail out on error)

%% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)

% Calculate time-domain pairwise-conditional causalities - this just requires
% the autocovariance sequence.

ptic('*** autocov_to_pwcgc... ');
F = autocov_to_pwcgc(G);
ptoc;

% Check for failed GC calculation

assert(~isbad(F,false),'GC calculation failed');

% Significance test using theoretical null distribution, adjusting for multiple
% hypotheses.

pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
sig  = significance(pval,alpha,mhtc);

% Plot time-domain causal graph, p-values and significance.

figure(2); clf;
subplot(1,3,1);
plot_pw(F);
title('Pairwise-conditional GC');
subplot(1,3,2);
plot_pw(pval);
title('p-values');
subplot(1,3,3);
plot_pw(sig);
title(['Significant at p = ' num2str(alpha)])

% For good measure we calculate Seth's causal density (cd) measure - the mean
% pairwise-conditional causality. We don't have a theoretical sampling
% distribution for this.

cd = mean(F(~isnan(F)));
fprintf('\ncausal density = %f\n',cd);

%% Granger causality calculation: frequency domain  (<mvgc_schema.html#3 |A14|>)

% Calculate spectral pairwise-conditional causalities at given frequency
% resolution - again, this only requires the autocovariance sequence.

ptic('\n*** autocov_to_spwcgc... ');
f = autocov_to_spwcgc(G,fres);
ptoc;

% Check for failed spectral GC calculation

assert(~isbad(f,false),'spectral GC calculation failed');

% Plot spectral causal graph.

figure(3); clf;
plot_spw(f,fs);

%% Granger causality calculation: frequency domain -> time-domain  (<mvgc_schema.html#3 |A15|>)

% Check that spectral causalities average (integrate) to time-domain
% causalities, as they should according to theory.

fprintf('\nchecking that frequency-domain GC integrates to time-domain GC... \n');
Fint = smvgc_to_mvgc(f); % integrate spectral MVGCs
mad = maxabs(F-Fint);
madthreshold = 1e-5;
if mad < madthreshold
    fprintf('maximum absolute difference OK: = %.2e (< %.2e)\n',mad,madthreshold);
else
    fprintf(2,'WARNING: high maximum absolute difference = %e.2 (> %.2e)\n',mad,madthreshold);
end
