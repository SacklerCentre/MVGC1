
subsample_parms;

% Get VAR coefficients for 5-node test network.

AT = var_specrad(var5_test,rho);
[nvars,~,morder] = size(AT);  % number of variables

if isempty(bsize), bsize = morder; end
nobs = floor(nobs/bsize)*bsize; % ensure permutation block size divides number of observations exactly

% Residuals covariance matrix.

rstate = rng_seed(cseed);
VT = random_covmat(nvars,grho);
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
fN = empirical_var_to_spwcgc(AT,VT,nobs,ntrials,fres,true,nsamps,'OLS');
rng_restore(rstate);

rstate = rng_seed(pseed);
fP = permtest_tsdata_to_spwcgc(X,morder,fres,bsize,nsamps,'OLS');
rng_restore(rstate);

if length(fband) == 1
    i = 1+round(fres*fband);
    FN = fN(:,:,:,i);
    FP = fP(:,:,:,i);
else
    FN = smvgc_to_mvgc(fN,fband);
    FP = smvgc_to_mvgc(fP,fband);
end

figure(1); clf
k = 0;
for i = 1:nvars
    for j = 1:nvars
        k = k+1;
        if i ~= j
            subplot(nvars,nvars,k);
            FNij = FN(:,i,j);
            FPij = FP(:,i,j);
            if isempty(Fmin), Fmn = min(min(FNij),min(FPij)); else Fmn = Fmin; end;
            if isempty(Fmax), Fmx = max(max(FNij),max(FPij)); else Fmx = Fmax; end;
            FVij = linspace(Fmn,Fmx,Fres)';
            PEij = empirical_cdf(FVij,FNij);
            PBij = empirical_cdf(FVij,FPij);
            PTij = mvgc_cdf(FVij,0,morder,nobs,ntrials,1,1,nvars-2,tstat);
            plot(FVij,[PTij PEij PBij]);
        end
    end
end
