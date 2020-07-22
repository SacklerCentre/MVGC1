%% MVGC toolbox demo
%
% Demonstrates typical usage of the MVGC routines on generated VAR data for
% a 5-node network with known causal structure (see image file
% 'var5_test.png'). Estimates a VAR model, outputs model statistcs and
% calculates time and frequency domain pairwise-conditional Granger
% causalities.

%% Parameters

% Note: you can set the regression mode ('OLS' or 'LWR') via 'set_regmode'. If
% not set explicitly, the default (OLS) is used.

ntrials   = 5;      % number of trials
nobs      = 1000;   % number of observations per trial

alpha     = 0.05;   % significance level for all statistical tests
stlags    = [];     % number of lags for unit-root stationarity test

%seed     = 81761;  % random seed (0 for unseeded)
seed      = 0;      % random seed (0 for unseeded)

%% Generate VAR data

% Seed random number generator.

rng_seed(seed);

% Get VAR coefficients for 5-node test network.

AT = var5_test;
nvars = size(AT,1); % number of variables

% Residuals covariance matrix.

SIGT = eye(nvars);

% Generate VAR time series data with normally distributed residuals for
% specified coefficients and covariance matrix.

ptic('\n*** var_to_tsdata... ');
X = var_to_tsdata(AT,SIGT,nobs,ntrials);
ptoc;

%% check stationarity

% Augmented Dickey-Fuller unit-root test

[adftstat,adfcval] = mvgc_adf(X,alpha,stlags);
fprintf('\nADF statistics (critical value = %f)\n',adfcval); disp(adftstat);
adfsig = adftstat > adfcval; % unit root; but how do we correct for multiple hypotheses?
adfnonstat = find(adfsig);
if isempty(adfnonstat)
    fprintf('all time series are stationary by ADF test at significance %g\n',alpha);
else
    if ntrials > 1 % multitrial
        for r = 1:ntrials
            adfnonstat = find(adfsig(r,:));
            if ~isempty(adfnonstat)
                fprintf(2,'WARNING: non-stationary time series by ADF test at significance %g for trial %d, variable(s): %s\n',alpha,r,num2str(adfnonstat));
            end
        end
    else
        fprintf(2,'WARNING: non-stationary time series by ADF test at significance %g for variable(s): %s\n',alpha,num2str(adfnonstat));
    end
end

% KPSS unit-root test

[ksstat,kscval] = mvgc_kpss(X,alpha,stlags);
fprintf('\nKPSS statistics (critical value = %f)\n',kscval); disp(ksstat);
kssig = ksstat > kscval; % unit root; but how do we correct for multiple hypotheses?
ksnonstat = find(kssig);
if isempty(ksnonstat)
    fprintf('all time series are stationary by KPSS test at significance %g\n',alpha);
else
    if ntrials > 1 % multitrial
        for r = 1:ntrials
            ksnonstat = find(kssig(r,:));
            if ~isempty(ksnonstat)
                fprintf(2,'WARNING: non-stationary time series by KPSS test at significance %g for trial %d, variable(s): %s\n',alpha,r,num2str(ksnonstat));
            end
        end
    else
        fprintf(2,'WARNING: non-stationary time series by KPSS test at significance %g for variable(s): %s\n',alpha,num2str(ksnonstat));
    end
end
