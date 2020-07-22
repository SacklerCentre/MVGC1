%% MVGC confidence intervals test demo
%
% Demonstrates calculation of confidence intervals with MVGC on generated VAR
% data for a 5-node network with known causal structure (see <var5_test.html
% |var5_test|>). Pairwise-conditional Granger causalities are estimated and
% confidence intervals constructed using both theoretical and empirical
% distributions.
%
% *_Note:_* the empirical MVGC distributions are calculated from data
% independently generated from the same model as used to calculate the test
% MVGC; in practice they might be aquired by bootstrap sampling (see
% <bootstrap_tsdata_to_mvgc.html |bootstrap_tsdata_to_mvgc|>).
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%% Parameters

ntrials   = 10;     % number of trials
nobs      = 100;    % number of observations per trial
nsamps    = 100;    % number of empirical samples

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)

acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)

tstat     = 'chi2'; % statistical test for MVGC: 'F' for Granger's F-test, 'chi2' for Geweke's chi2 test or leave empty for default
alpha     = 0.05;   % significance level for all statistical tests
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

seed      = 0;      % random seed (0 for unseeded)

%% Simulation

% Seed random number generator.

rng_seed(seed);

% Get VAR coefficients for 5-node test network.

AT = var5_test;
nvars = size(AT,1);  % number of variables

% Residuals covariance matrix.

SIGT = eye(nvars);

morder = size(AT,3); % actual model order; on real data - i.e. with no generative model
                     % available - use information criteria to estimate (see 'mvgc_demo')

ns1 = nsamps+1; % first sample will be MVGC to test
F = zeros(ns1,nvars,nvars);
s = 1;
r = 0;
while s <= nsamps+1
    fprintf('sample %d of %d',s,ns1);

    % Generate VAR time series data with normally distributed residuals for
    % specified coefficients and covariance matrix.

    X = var_to_tsdata(AT,SIGT,nobs,ntrials);

    % Calculate VAR model

    [A,SIG] = tsdata_to_var(X,morder,regmode);

    % Check for failed regression.

    if isbad(A)
        fprintf(' *** VAR estimation failed\n');
        r = r+1;
    else

        % Calculate autocovariance.

        [G,info] = var_to_autocov(A,SIG,acmaxlags);

        % Check for errors.

        if info.error
            fprintf(' *** bad VAR: %s\n',info.errmsg);
            r = r+1;
        else

            % Calculate time-domain pairwise conditional causalities.

            F(s,:,:) = autocov_to_pwcgc(G);
            if isbad(F(s,:,:),false)
                fprintf(2,' *** GC calculation failed\n');
                r = r+1;
            else
                fprintf('\n');
                s = s+1;
            end

        end
    end
end
fprintf('re-tries = %d\n',r);

FEMP = F(2:end,:,:);   % MVGC empirical distribution
F = squeeze(F(1,:,:)); % MVGC to test

%% Confidence intervals

% Estimate confidence intervals for F from theoretical and empirical distributions.

[FTUP,FTLO] = mvgc_confint(alpha,F,morder,nobs,ntrials,1,1,nvars-2,tstat);

[FEUP,FELO] = empirical_confint(alpha,FEMP);

% Critical GC value - theoretical only, as we haven't calculated an empirical
% null distribution.

FTCRIT = mvgc_cval(alpha,morder,nobs,ntrials,1,1,nvars-2,tstat);

% Plot them

figure(1); clf
subplot(2,1,1);
plot_confints(F,FTUP,FTLO,FTCRIT);
title(sprintf('Theoretical distribution\nconfidence intervals at alpha = %g',alpha));
subplot(2,1,2);
plot_confints(F,FEUP,FELO);
title(sprintf('Empirical distribution\nconfidence intervals at alpha = %g',alpha));

%%
% <mvgc_demo_confint.html back to top>
