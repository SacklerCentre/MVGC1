%% MVGC demo: state-space method.
%
% This demo demonstrates basic usage of the MVGC toolbox using the state-space
% method introduced in [4]. We recommend this approach as default, since it is
% in general faster and potentially more accurate. For the current default
% autocovariance-based method, see <mvgc_demo.html |mvgc_demo|>.
%
% Demonstrates typical usage of the MVGC toolbox on generated VAR data for a
% 5-node network with known causal structure (see <var5_test.html |var5_test|>).
% Estimates a VAR model and calculates time- and frequency-domain
% pairwise-conditional Granger causalities (also known as the "causal graph").
%
% This script is a good starting point for learning the MVGC approach to
% Granger-causal estimation and statistical inference. It may serve as a useful
% template for your own code.
%
% *_FAQ:_* _Why do the spectral causalities look so smooth?_
%
% This is because spectral quantities are calculated from the estimated VAR,
% rather than sampled directly. This is in accordance with the MVGC design
% principle that all causal estimates be based on the <mvgc_demo.html#6
% estimated VAR model> for your data, and guarantees that spectral causalities
% <mvgc_demo.html#10 integrate correctly> to time-domain causality as theory
% requires. See [1] for details.
%
% *_Note_*: Do _not_ pre-filter your data prior to GC estimation, _except_
% possibly to improve stationarity (e.g notch-filtering to eliminate line noise
% or high-pass filtering to suppress low-frequency transients). Pre-filtering
% (of stationary data) may seriously degrade Granger-causal inference! If you
% want (time-domain) GC over a limited frequency range, rather calculate
% "band-limited" GC; to do this, calculate frequency-domain GCs over the full
% frequency range, then integrate over the desired frequency band [3]; see
% <smvgc_to_mvgc.html |smvgc_to_mvgc|>.
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% [2] A. B. Barrett, L. Barnett and A. K. Seth, "Multivariate Granger causality
% and generalized variance", _Phys. Rev. E_ 81(4), 2010.
%
% [3] L. Barnett and A. K. Seth, "Behaviour of Granger causality under
% filtering: Theoretical invariance and practical application", _J. Neurosci.
% Methods_ 201(2), 2011.
%
% [4] L. Barnett and A. K. Seth, "Granger causality for state-space models",
% _Phys. Rev. E 91(4) Rapid Communication_, 2015
% [ <matlab:open('ssgc_preprint.pdf') preprint> ].
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%% Parameters

ntrials   = 10;     % number of trials
nobs      = 1000;   % number of observations per trial

regmode   = 'LWR';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 20;     % maximum model order for model order estimation

acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)

tstat     = 'F';    % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDRD'; % multiple hypothesis test correction (see routine 'significance')

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

%% Model order estimation (<mvgc_schema.html#3 |A2|>)

% Calculate information criteria up to specified maximum model order.

ptic('\n*** tsdata_to_infocrit\n');
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
ptoc('*** tsdata_to_infocrit took ');

% Plot information criteria.

figure(1); clf;
plot_tsdata([AIC BIC]',{'AIC','BIC'},1/fs);
title('Model order estimation');

amo = size(AT,3); % actual model order

fprintf('\nbest model order (AIC) = %d\n',moAIC);
fprintf('best model order (BIC) = %d\n',moBIC);
fprintf('actual model order     = %d\n',amo);

% Select model order.

if     strcmpi(morder,'actual')
    morder = amo;
    fprintf('\nusing actual model order = %d\n',morder);
elseif strcmpi(morder,'AIC')
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
assert(~isbad(A),'VAR estimation failed - bailing out');
ptoc;

% Report information on the estimated VAR, and check for errors.
%
% _IMPORTANT:_ We check the VAR model for stability and symmetric
% positive-definite residuals covariance matrix. _THIS CHECK SHOULD ALWAYS BE
% PERFORMED!_ - subsequent routines may fail if there are errors here. If there
% are problems with the data (e.g. non-stationarity, colinearity, etc.) there's
% also a good chance they'll show up at this point - and the diagnostics may
% supply useful information as to what went wrong.

info = var_info(A,SIG);
assert(~info.error,'VAR error(s) found - bailing out');

%% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)

% Calculate time-domain pairwise-conditional causalities from VAR model parameters
% by state-space method [4]. The VAR model is transformed into an equivalent state-
% space model for computation. Also return p-values for specified test (F-test or
% likelihood-ratio test; this is optional - if p-values are not required, then it
% is not necessary to supply time series |X|, regression mode |regmode|, or test
% specification |tstat|).

ptic('*** var_to_pwcgc... ');
[F,pval] = var_to_pwcgc(A,SIG,X,regmode,tstat);
ptoc;

% Check for failed GC calculation

assert(~isbad(F,false),'GC calculation failed - bailing out');

% Significance-test p-values, correcting for multiple hypotheses.

sig = significance(pval,alpha,mhtc);

% Plot time-domain causal graph, p-values and significance.

figure(2); clf;
sgtitlex('Pairwise-conditional Granger causality - time domain');
subplot(1,3,1);
plot_pw(F);
title('Pairwise-conditional GC');
subplot(1,3,2);
plot_pw(pval);
title(['p-values (' tstat '-test)']);
subplot(1,3,3);
plot_pw(sig);
title(['Significant at \alpha = ' num2str(alpha)]);

%% Granger causality calculation: frequency domain  (<mvgc_schema.html#3 |A14|>)

% If not specified, we set the frequency resolution to something sensible. Warn if
% resolution is very large, as this may lead to excessively long computation times,
% and/or out-of-memory issues.

if isempty(fres)
    fres = 2^nextpow2(info.acdec); % based on autocorrelation decay; alternatively, you could try fres = 2^nextpow2(nobs);
	fprintf('\nfrequency resolution auto-calculated as %d (increments ~ %.2gHz)\n',fres,fs/2/fres);
end
if fres > 20000 % adjust to taste
	fprintf(2,'\nWARNING: large frequency resolution = %d - may cause computation time/memory usage problems\nAre you sure you wish to continue [y/n]? ',fres);
	istr = input(' ','s'); if isempty(istr) || ~strcmpi(istr,'y'); fprintf(2,'Aborting...\n'); return; end
end

% Calculate spectral pairwise-conditional causalities at given frequency
% resolution by state-space method.

ptic('\n*** var_to_spwcgc... ');
f = var_to_spwcgc(A,SIG,fres);
assert(~isbad(f,false),'spectral GC calculation failed - bailing out');
ptoc;

% Plot spectral causal graph.

figure(3); clf;
sgtitlex('Pairwise-conditional Granger causality - frequency domain');
plot_spw(f,fs);

%% Granger causality calculation: frequency domain -> time-domain  (<mvgc_schema.html#3 |A15|>)

% Check that spectral causalities average (integrate) to time-domain
% causalities, as they should according to theory.

fprintf('\nfrequency-domain GC integration check... ');
Fint = smvgc_to_mvgc(f); % integrate spectral MVGCs
amax = maxabs(F+Fint)/2;
if amax < 1e-5; amax = 1; end % in case all GCs very small
mre = maxabs(F-Fint)/amax;
if mre < 1e-5
    fprintf('OK (maximum relative error ~ %.0e)\n',mre);
else
    fprintf(2,'WARNING: high maximum relative error ~ %.0e\n',mre);
end

%%
% <mvgc_demo_statespace.html back to top>
