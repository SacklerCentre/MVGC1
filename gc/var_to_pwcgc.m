function F = var_to_pwcgc(A,SIG)

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(SIG);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

DV = diag(SIG);
LDSIG = log(DV);

F = nan(n);
for j = 1:n
    jo = [1:j-1 j+1:n]; % omit j
	[~,SIGj,rep] = var2riss(A,SIG,j,jo); % residuals covariance matrix of reduced model
    if sserror(rep,j), continue; end     % check DARE report, bail out on error
    F(jo,j) = log(diag(SIGj))-LDSIG(jo);
end
