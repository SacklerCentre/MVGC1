function F = mvgc_cval_old(alpha,morder,nobs,ntrials,nxvars,nyvars,nzvars,tstat)

% F = mvgc_cval(alpha,morder,nobs,ntrials,nxvars,nyvars,nzvars,tstat)
%
% Return critical MVGC value at significance alpha for time-domain sample MVGC,
% based on theoretical (asymptotic) null distributions.
%
% morder is model order (number of lags), ntrials is number of trials, nobs is
% number of observations per trial, nxvars is number of predictee variables,
% nyvars is number of predictor variables, nzvars is number of conditioning
% variables. NaNs are ignored.
%
% The statistical test is specified by the tstat parameter, which may be 'F'
% (for Granger's F-test) or 'chi2' (for Geweke's chi2 test); however, for a
% multivariate predictee (nxvars > 1) only the chi2 test is available. Both
% tests are asymptotic, although the F-test (where available) may be more
% accurate for small samples; hence the default is for the F-test if nxvars == 1
% and the chi2 test otherwise.
%
% IMPORTANT: To test for significance you should adjust alpha for multiple null
% hypotheses; see routines 'significance' and also 'mvgc_pval'.

if nargin < 7 || isempty(nzvars), nzvars = 0; end % unconditional

if nargin < 8 || isempty(tstat);
    ftest = nxvars == 1;                  % default: use F-test for univariate predictee, chi2 for multivariate
else
    switch lower(tstat)
        case 'f'    % Granger F-test
            assert(nxvars == 1,'Cannot use F-test for multivariate predictee');
            ftest = true;
        case 'chi2' % Geweke chi2 test
            ftest = false;
        otherwise
            error('unknown statistical test');
    end
end

enobs = ntrials*(nobs-morder);            % effective number of observations: note that a p-lag
                                          % autoregression loses morder observations per trial
if ftest % Granger F-test

    df1 = morder*nyvars;                  % #{full model parameters} - #{reduced model parameters}
    df2 = enobs-morder*(1+nyvars+nzvars); % #{observations} - #{full model parameters}
    F = log(1+(df1/df2)*finv(1-alpha,df1,df2));

else     % Geweke chi2 test

    df = morder*nxvars*nyvars;            % note that df does not depend on the number of conditioning variables
    F = chi2inv(1-alpha,df)/enobs;

end
