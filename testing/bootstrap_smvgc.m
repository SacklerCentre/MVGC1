
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
fE = empirical_var_to_smvgc(AT,VT,nobs,ntrials,x,y,fres,false,nsamps,'OLS');
rng_restore(rstate);

rstate = rng_seed(bseed);
fB = bootstrap_tsdata_to_smvgc(X,x,y,morder,fres,nsamps);
rng_restore(rstate);

if length(fband) == 1
    i = 1+round(fres*fband);
    FE = fE(:,i);
    FB = fB(:,i);
else
    FE = smvgc_to_mvgc(fE,fband);
    FB = smvgc_to_mvgc(fB,fband);
end

if isempty(Fmin), Fmin = min(min(FE),min(FB)); end;
if isempty(Fmax), Fmax = max(max(FE),max(FB)); end;
FVALS = linspace(Fmin,Fmax,Fres)';

PE = empirical_cdf(FVALS,FE);
PB = empirical_cdf(FVALS,FB);

F = nanmean(FE);

PT = mvgc_cdf(FVALS,F,morder,nobs,ntrials,nx,ny,nz,tstat);

plot(FVALS,[PT PE PB]);
legend('theoretical','empirical','bootstrap','Location','NorthEastOutside');
