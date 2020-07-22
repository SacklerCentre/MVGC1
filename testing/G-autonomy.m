%% G-autonomy.m
%
% This script demonstrating G-autonomy calculation is adapted from the main MVGC
% demo script, mvgc_demo.m
%
%% Parameters

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

morder    = 'AIC';  % model order to use ('AIC', 'BIC' or supplied numerical value)
momax     = 20;     % maximum model order for model order estimation

acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)

tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

fs        = 250;    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)

%% Load your data

X = ???; % the (multi-)variable X you want to test for autonomy (can be multi-trial)
Z = ???; % the external variables

U = cat(1,X,Z); % the "universe" of all variables

[nvars,nobs,ntrials] = size(U);
nxvars = size(X,1); % dimension of X
x = 1:nxvars;       % indices of X in U


%% Model order estimation

% Calculate information criteria up to specified maximum model order.

ptic('\n*** tsdata_to_infocrit\n');
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(U,momax,icregmode);
ptoc('*** tsdata_to_infocrit took ');

% Plot information criteria.

figure(1); clf;
plot_tsdata([AIC BIC]',{'AIC','BIC'},1/fs);
title('Model order estimation');

fprintf('\nbest model order (AIC) = %d\n',moAIC);
fprintf('best model order (BIC) = %d\n',moBIC);

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

%% G-autonomy calculation: time domain  (<mvgc_schema.html#3 |A13|>)

% Full regression

ptic('\n*** tsdata_to_var... ');
[~,SIG] = tsdata_to_var(U,morder,regmode);
ptoc;

% Check for failed regression

assert(~isbad(A),'VAR estimation failed');

% Calculate time-domain pairwise-conditional causalities - this just requires
% the autocovariance sequence.

% full regression

[A,SIG,E] = tsdata_to_var(U,p,regmode)

tsdata_to_var
owstate = warn_supp;
[~,SIG ] = autocov_to_var(G);
warn_test(owstate,    'in full regression - bad autocovariance matrix? Check output of ''var_info''');
if warn_if(isbad(SIG),'in full regression - regression failed'), return; end % show-stopper!

% reduced regression

owstate = warn_supp;
[~,SIGR] = autocov_to_var(G(xz,xz,:));    % reduced regression
warn_test(owstate,     'in reduced regression - bad autocovariance matrix? Check output of ''var_info''');
if warn_if(isbad(SIGR),'in reduced regression - regression failed'), return; end % show-stopper!

x = 1:length(x);
F = log(det(SIGR(x,x)))-log(det(SIG(x,x)));









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

%%
% <mvgc_demo.html back to top>
