%% consistency
%
% Calculate VAR model consistency
%
% <matlab:open('consistency.m') code>
%
%% Syntax
%
%     cons = consistency(X,E)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     X          multi-trial time series data
%     E          residuals time series
%
% _output_
%
%     cons       consistency statistic
%
%% Description
%
% Check what proportion of the correlation structure in data is accounted
% for by a VAR model. The data |X| and residuals |E| may be single- or
% multi-trial time series. A value of |cons > 0.8| might be considered
% to indicate reasonable consistency.
%
%% References
%
% [1] M. Ding, S. L. Bressler, W. Yang and H. Liang, "Short-window spectral
% analysis of cortical event-related potentials by adaptive multivariate
% autoregressive modeling: data preprocessing, model validation, and
% variability assessment", _Biol. Cybern._, 83, 2000.
%
% [2] A. K. Seth, "A MATLAB toolbox for Granger causal connectivity
% analysis", _Journal of Neuroscience Methods_ 186, 2010.
%
%% See also
%
% <tsdata_to_var.html |tsdata_to_var|> |
% <mvgc_demo_stats.html |mvgc_demo_stats|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function cons = consistency(X,E)

[n,m,N] = size(X);
p = m-size(E,2);     % number of lags in model
assert(m >= n,'too few observations');
assert(size(E,1) == n && size(E,3) == N,'residuals don''t match data');
assert(p > 0,'bad number of lags');

X = demean(X);

M = N*(m-p);                     % effective number of observations
X = reshape(X(:,p+1:m,:),n,M);   % concatenate trials for data
E = E(:,:);                      % concatenate trials for residuals

Y = X-E;                         % prediction

Rr = (X*X')/(M-1);               % covariance estimate
Rs = (Y*Y')/(M-1);               % covariance estimate

cons = 1 - norm(Rs-Rr)/norm(Rr); % compare matrix norms
