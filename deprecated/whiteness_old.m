function [dw,pval] = whiteness_old(X,E)

% [dw,pval] = whiteness(X,E)
%
% Check for uncorrelated residuals using Durbin-Watson test. X is the data,
% E the residuals. The data X and residuals E may be single- or multi-trial
% time series. The Durbin-Watson test statistic for each variable is
% returned in dw, with p-values in pval.
%
% NOTE: To test for significance you should correct for multiple null
% hypotheses, or test for false discovery rate; see routine significance.m.
%
% TODO: This is not really satisfactory for multi-trial data. We should
% investigate:
%
% Bhargava, Alok, Franzini, L., Narendranathan, W. (1982): "Serial
% Correlation and the Fixed Effects Model". Review of Economic Studies, 49,
% p. 533?549.
%
% A permutation test for whiteness would be nice too...

[n,m,N] = size(X);

assert(m >= n,'too few observations');
assert(size(E,1) == n,'residuals don''t match data');
assert(size(E,2) < m,'bad number of lags');

X = demean(X);

dw = zeros(N,n);
pval = zeros(N,n);
for r = 1:N % note: this works for single-trial (N == 1) too
    for i = 1:n
        [pval(r,i),dw(r,i)] = dwtest(E(i,:,r)',X(:,:,r)'); % Evaluate DW statistics for each residuals series
    end
end
