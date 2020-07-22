
subsample_parms;

assert(causal,'must be causal!');

% Get VAR coefficients for 9-node test network.

AT = var_specrad(var9_test,rho);
[nvars,~,morder] = size(AT);  % number of variables

if isempty(bsize), bsize = morder; end
nobs = floor(nobs/bsize)*bsize; % ensure permutation block size divides number of observations exactly

nx = length(x);
ny = length(y);
nz = nvars-nx-ny;

% Residuals covariance matrix.

rstate = rng_seed(corrseed);
VT = random_covmat(nvars,corrfac);
rng_restore(rstate);
assert(isposdef(VT));
                     
rstate = rng_seed(xseed);
X = var_to_tsdata(AT,VT,nobs,ntrials);
rng_restore(rstate);

if isempty(fres);
    A = tsdata_to_var(X,morder,'OLS');
    rho = var_specrad(A);
    assert(rho < 1,'unstable VAR (unit root)');
    fres = aclags(rho)+1;
    fprintf('fres = %d\n\n',fres);
end

rstate = rng_seed(eseed);
fN = empirical_var_to_smvgc(AT,VT,nobs,ntrials,x,y,fres,true,nsamps,'OLS');
rng_restore(rstate);

rstate = rng_seed(pseed);
fP = permtest_tsdata_to_smvgc(X,x,y,morder,fres,bsize,nsamps,'OLS');
rng_restore(rstate);

if length(fband) == 1
    i = 1+round(fres*fband);
    FN = fN(:,i);
    FP = fP(:,i);
else
    FN = smvgc_to_mvgc(fN,fband);
    FP = smvgc_to_mvgc(fP,fband);
end

if isempty(Fmin), Fmin = min(min(FN),min(FP)); end;
if isempty(Fmax), Fmax = max(max(FN),max(FP)); end;
FVALS = linspace(Fmin,Fmax,Fres)';

PN = empirical_cdf(FVALS,FN);
PP = empirical_cdf(FVALS,FP);

PT = mvgc_cdf(FVALS,0,morder,nobs,ntrials,nx,ny,nz,tstat);

plot(FVALS,[PT PN PP]);
legend('theoretical','null','permutation','Location','NorthEastOutside');
