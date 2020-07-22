
subsample_parms;

% Get VAR coefficients for 5-node test network.

AT = var_specrad(var5_test,rho);
[nvars,~,morder] = size(AT);  % number of variables

% Residuals covariance matrix.

rstate = rng_seed(cseed);
VT = random_covmat(nvars,grho);
rng_restore(rstate);
assert(isposdef(VT));

rstate = rng_seed(eseed);
FE = empirical_var_to_pwcgc(AT,VT,nobs,ntrials,false,nsamps,'OLS');
rng_restore(rstate);
                     
rstate = rng_seed(xseed);
X = var_to_tsdata(AT,VT,nobs,ntrials);
rng_restore(rstate);

rstate = rng_seed(bseed);
FB = bootstrap_tsdata_to_pwcgc(X,morder,nsamps);
rng_restore(rstate);

figure(1); clf
k = 0;
for i = 1:nvars
    for j = 1:nvars
        k = k+1;
        if i ~= j
            subplot(nvars,nvars,k);
            FEij = FE(:,i,j);
            FBij = FB(:,i,j);
            Fij = nanmean(FEij);
            if isempty(Fmin), Fmn = min(min(FEij),min(FBij)); else Fmn = Fmin; end;
            if isempty(Fmax), Fmx = max(max(FEij),max(FBij)); else Fmx = Fmax; end;
            FVij = linspace(Fmn,Fmx,Fres)';
            PEij = empirical_cdf(FVij,FEij);
            PBij = empirical_cdf(FVij,FBij);
            PTij = mvgc_cdf(FVij,Fij,morder,nobs,ntrials,1,1,nvars-2,tstat);
            plot(FVij,[PTij PEij PBij]);
        end
    end
end;
