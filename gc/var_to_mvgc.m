%% var_to_mvgc
%
% Calculate conditional time-domain MVGC (multivariate Granger causality)
%
% <matlab:open('var_to_mvgc.m') code>
%
%% Syntax
%
%     [F,pval] = var_to_mvgc(A,SIG,x,y,X,regmode,tstat)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     A          VAR coefficients matrix
%     SIG        residuals covariance matrix
%     x          vector of indices of target (causee) multi-variable
%     y          vector of indices of source (causal) multi-variable
%     X          multi-trial time series data
%     regmode    regression mode: 'LWR' or 'OLS'
%     tstat      statistical inference test: 'F' for F-test, or 'LR' for likelihood-ratio (chi^2) test
%
% _output_
%
%     F          Granger causality
%     pval       p-value for specified statistical test
%
%% Description
%
% Returns the time-domain MVGC
%
% <<eq_mvgc.png>>
%
% from the variable |Y| (specified by the vector of indices |y|) to the
% variable |X| (specified by the vector of indices |x|), conditional on all
% other variables |Z| represented in |A| and |SIG|, for a stationary VAR process
% with VAR coefficients matrix |A| and residuals covariance matrix |SIG| - see
% ref. [1].
%
% The algorithm first converts the VAR parameters to state-space innovations form
% (see <var2riss.html |var2riss|>) then applies the method detailed in ref. [2]
% to calculate Granger causality |F| from the state-space parameters.
%
% If the return value |pval| is specified, then p-values are calculated for the
% null hypothesis of zero Granger causality according to an F- or chi^2 test. In
% this case, the parameters |X|, |regmode| and |tstat| must be supplied.
%
% The caller should take note of any warnings issued by this function and test
% results with a call <isbad.html |isbad|> (|F|, |false|).
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% [2] L. Barnett and A. K. Seth, "Granger causality for state-space models",
% _Phys. Rev. E 91(4) Rapid Communication_, 2015
% [ <matlab:open('ssgc_preprint.pdf') preprint> ].
%
%% See also
%
% <var2riss.html |var2riss|> |
% <isbad.html |isbad|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%
function [F,pval] = var_to_mvgc(A,SIG,x,y,X,regmode,tstat)

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(SIG);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match VAR coefficients matrix');

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(length(unique([x y])) == length([x y]),'x and y indices must be unique and non-overlapping');
assert(all(x >=1 & x <= n),'Some x indices out of range');
assert(all(y >=1 & y <= n),'Some y indices out of range');

z  = 1:n; z([x y]) = []; % indices of other variables (to condition out)
xz = [x z];              % indices of variables in reduced model (omit source variables)

nx = length(x);
ny = length(y);
xr = 1:nx; % indices of x in reduced model

F = NaN;
[~,SIGR,rep] = var2riss(A,SIG,y,xz); % residuals covariance matrix of reduced model
if sserror(rep), return; end         % check DARE report, bail out on error
F = logdet(SIGR(xr,xr)) - logdet(SIG(x,x));

if nargout > 1 % calculate stats
	assert(nargin > 5, 'Must supply regression mode for stats (same mode as used for VAR model estimate)');
	assert(~isempty(X),'Must supply time-series data for stats');
	if nargin < 7 || isempty(tstat), tstat = 'F'; end % default is F-test (better for shorter time series)
	[n1,m,N] = size(X);
	assert(n1 == n,    'Time series does not match VAR coefficients matrix');
	M  = N*(m-p); % effective number of observations
	d  = p*nx*ny; % degrees of freedom
	[~,VR]  = tsdata_to_var(X(xz,:,:),p,regmode);  % reduced regression
	if     strcmpi(tstat,'F')  % F-test
		d2 = nx*(M-p*n)-1; % F df2
		K  = d2/d;         % F scaling factor
		stat  = trace(SIGR(xr,xr))/trace(SIG(x,x)) - 1; % F-test statistic
		pval = 1-fcdf(K*stat,d,d2);
	elseif strcmpi(tstat,'LR') % Likelihood-ratio test
		stat = logdet(SIGR(xr,xr)) - logdet(SIG(x,x)); % likelihood-ratio test statistic
		pval = 1-chi2cdf(M*stat,d);
	else
		error('Unknown statistical test');
	end
end
