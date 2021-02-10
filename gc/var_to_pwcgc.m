function [F,pval] = var_to_pwcgc(A,SIG,X,regmode,tstat)

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(SIG);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

DSIG = diag(SIG);
LDSIG = log(DSIG);

F = nan(n);
for j = 1:n
    jo = [1:j-1 j+1:n]; % omit j
	[~,SIGj,rep] = var2riss(A,SIG,j,jo); % residuals covariance matrix of reduced model
    if sserror(rep,j), continue; end     % check DARE report, bail out on error
    F(jo,j) = log(diag(SIGj))-LDSIG(jo);
end

if nargout > 1 % calculate stats
	assert(nargin > 3, 'Must supply regression mode for stats (same mode as used for VAR model estimate)');
	assert(~isempty(X),'Must supply time-series data for stats');
	if nargin < 5 || isempty(tstat), tstat = 'F'; end % default is F-test (better for shorter time series)
	[n1,m,N] = size(X);
	assert(n1 == n,    'Time series does not match VAR coefficients matrix');
    M  = N*(m-p);  % effective number of observations
    d  = p;        % degrees of freedom
	stat = nan(n);
	if     strcmpi(tstat,'F')  % F-test
		d2 = M-p*n-1;  % F df2
		K  = d2/d;     % F scaling factor
		for j = 1:n
			jo = [1:j-1 j+1:n]; % omit j
			[~,SIGj] = tsdata_to_var(X(jo,:,:),p,regmode); % reduced regression
			DSIGj = diag(SIGj);
			stat(jo,j) = diag(SIGj)./DSIG(jo) - 1; % F-test statistic
		end
		pval = 1-fcdf(K*stat,d,d2);
	elseif strcmpi(tstat,'LR') % Likelihood-ratio test
		for j = 1:n
			jo = [1:j-1 j+1:n]; % omit j
			[~,SIGj] = tsdata_to_var(X(jo,:,:),p,regmode); % reduced regression
			stat(jo,j) = log(diag(SIGj)) - LDSIG(jo);     % likelihood-ratio test statistic
		end
		pval = 1-chi2cdf(M*stat,d);
	else
		error('Unknown statistical test');
	end
end
