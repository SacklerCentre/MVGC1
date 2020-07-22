
subsample_parms;

% Get VAR coefficients for 9-node test network.

AT = var_specrad(var9_test,rho);
[nvars,~,morder] = size(AT);  % number of variables

% Residuals covariance matrix.

rstate = rng_seed(cseed);
VT = random_covmat(nvars,grho);
rng_restore(rstate);
assert(isposdef(VT));

nx = length(x);
ny = length(y);
nz = nvars-nx-ny;

rstate = rng_seed(eseed);
FE = empirical_var_to_mvgc(AT,VT,nobs,ntrials,x,y,false,nsamps,'OLS');
rng_restore(rstate);
                     
rstate = rng_seed(xseed);
U = var_to_tsdata(AT,VT,nobs,ntrials);
rng_restore(rstate);

rstate = rng_seed(bseed);
FB = bootstrap_tsdata_to_mvgc(U,x,y,morder,nsamps);
rng_restore(rstate);

if isempty(Fmin), Fmin = min(min(FE),min(FB)); end;
if isempty(Fmax), Fmax = max(max(FE),max(FB)); end;
FVALS = linspace(Fmin,Fmax,Fres)';

PE = empirical_cdf(FVALS,FE);
PB = empirical_cdf(FVALS,FB);

F = nanmean(FE);

PT = mvgc_cdf(FVALS,F,morder,nobs,ntrials,nx,ny,nz,tstat);

plot(FVALS,[PT PE PB]);
legend('theoretical','empirical','bootstrap','Location','NorthEastOutside');
