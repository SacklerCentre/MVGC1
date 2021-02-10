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
xz = [x z];               % indices of variables in reduced model (omit source variables)

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
	[~,VR]  = tsdata_to_var(X(r,:,:),p,regmode);  % reduced regression
	if     strcmpi(tstat,'F')  % F-test
		d2 = nx*(M-p*n)-1; % F df2
		K  = d2/d;         % F scaling factor
		stat  = trace(VR(xr,xr))/trace(V(x,x)) - 1; % F-test statistic
		pval = 1-fcdf(K*stat,d,d2);
	elseif strcmpi(tstat,'LR') % Likelihood-ratio test
		stat = logdet(VR(xr,xr)) - logdet(V(x,x)); % likelihood-ratio test statistic
		pval = 1-chi2cdf(M*stat,d);
	else
		error('Unknown statistical test');
	end
end
