
subsample_parms;

% Get VAR coefficients for 5-node test network.

AT = var_specrad(var5_test,rho);
[nvars,~,morder] = size(AT);  % number of variables

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
fE = empirical_var_to_spwcgc(AT,VT,nobs,ntrials,fres,false,nsamps,'OLS');
rng_restore(rstate);

rstate = rng_seed(bseed);
fB = bootstrap_tsdata_to_spwcgc(X,morder,fres,nsamps);
rng_restore(rstate);

if length(fband) == 1
    i = 1+round(fres*fband);
    FE = fE(:,:,:,i);
    FB = fB(:,:,:,i);
else
    FE = smvgc_to_mvgc(fE,fband);
    FB = smvgc_to_mvgc(fB,fband);
end

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
end
