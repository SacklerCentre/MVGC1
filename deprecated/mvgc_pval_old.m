function pval = mvgc_pval_old(F,morder,nobs,ntrials,nxvars,nyvars,nzvars,tstat)

% pval = mvgc_pval(F,morder,nobs,ntrials,nxvars,nyvars,nzvars,tstat)
%
% Return p-value for time-domain sample MVGC, based on theoretical
% (asymptotic) null distributions.
%
% F is sample MVGC (F may be a vector or matrix), morder is model order
% (number of lags), ntrials is number of trials, nobs is number of
% observations per trial, nxvars is number of predictee variables, nyvars
% is number of predictor variables, nzvars is number of conditioning
% variables. NaNs are ignored.
%
% The statistical test is specified by the tstat parameter, which may be
% 'F' (for Granger's F-test) or 'chi2' (for Geweke's chi2 test); however,
% for a multivariate predictee (nxvars > 1) only the chi2 test is available.
% Both tests are asymptotic, although the F-test (where available) may be
% more accurate for small samples; hence the default is for the F-test if
% nxvars == 1 and the chi2 test otherwise.
%
% IMPORTANT: To test for significance you should correct for multiple null
% hypotheses, or test for false discovery rate; see routine 'significance'.

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

pval = NaN(size(F));                      % output p-value matrix is same shape as F matrix
nn   = ~isnan(F);                         % indices of non-NaN F values (logical array)
F    = F(nn);                             % vectorise the non-NaN F values

enobs = ntrials*(nobs-morder);            % effective number of observations: note that a p-lag
                                          % autoregression loses morder observations per trial
if ftest % Granger F-test

    G = exp(F)-1;                         % the Granger form: (RSS_reduced - RSS_full) / RSS_full
    df1 = morder*nyvars;                  % #{full model parameters} - #{reduced model parameters}
    df2 = enobs-morder*(1+nyvars+nzvars); % #{observations} - #{full model parameters}
    pvalnn = 1-fcdf((df2/df1)*G,df1,df2);

else     % Geweke chi2 test

    df = morder*nxvars*nyvars;            % note that df does not depend on the number of conditioning variables
    pvalnn = 1-chi2cdf(enobs*F,df);

end

pval(nn) = pvalnn;                        % NaNs will be in same positions as they were in original F matrix
