
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

rstate = rng_seed(cseed);
VT = random_covmat(nvars,grho);
rng_restore(rstate);
assert(isposdef(VT));

rstate = rng_seed(eseed);
FN = empirical_var_to_mvgc(AT,VT,nobs,ntrials,x,y,true,nsamps,'OLS');
rng_restore(rstate);
                     
rstate = rng_seed(xseed);
U = var_to_tsdata(AT,VT,nobs,ntrials);
rng_restore(rstate);

rstate = rng_seed(pseed);
FP = permtest_tsdata_to_mvgc(U,x,y,morder,bsize,nsamps,'OLS');
rng_restore(rstate);

if isempty(Fmin), Fmin = min(min(FN),min(FP)); end;
if isempty(Fmax), Fmax = max(max(FN),max(FP)); end;
FVALS = linspace(Fmin,Fmax,Fres)';

PN = empirical_cdf(FVALS,FN);
PP = empirical_cdf(FVALS,FP);

PT = mvgc_cdf(FVALS,0,morder,nobs,ntrials,nx,ny,nz,tstat);

figure(1); clf
plot(FVALS,[PT PN PP]);
legend('theoretical','null','permutation','Location','NorthEastOutside');
