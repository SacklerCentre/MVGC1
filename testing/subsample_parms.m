ntrials   = 10;     % number of trials
nobs      = 1000;   % number of observations per trial
nsamps    = 100;    % number of empirical samples
bsize     = [];     % permutation block size (default to model order)

rho       = 0.8;    % spectral radius
grho      = 0.5;    % residuals approx generalised correlation
causal    = true;   % causal?

fs        = 200;    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)
fband     = [NaN 0.1];

tstat     = 'chi2'; % statistical test for MVGC: 'fN' for Granger's fN-test, 'chi2' for Geweke's chi2 test or leave empty for default
alpha     = 0.05;   % significance level for all statistical tests
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

cseed     = 8271;   % residuals covariance random seed (0 for unseeded)
eseed     = 6511;   % empirical random seed (0 for unseeded)
xseed     = 1782;   % time series random seed (0 for unseeded)
pseed     = 3192;   % bootsdtrap random seed (0 for unseeded)

Fmin = 0;
Fmax = [];
Fres = 200;

if causal
    x  = 5;
    y  = 3;
else
    x  = 3;
    y  = 5;
end
