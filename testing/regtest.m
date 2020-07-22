ntests    = 100;

ntrials   = 10;    % number of trials (if > 1 uses multi-trial routines)
nobs      = 1000;  % number of observations per trial

%seed      = 0;      % random seed (0 for unseeded)
seed      = 81761;  % random seed (0 for unseeded)

morder    = 10;
nvars     = 9;
rho       = 0.9;
corrfac   = 0.5;    % residuals correlation factor

rng_seed(seed);

dfaAO = zeros(ntests,1);
dfaAL = zeros(ntests,1);
dfaOL = zeros(ntests,1);

dfvAO = zeros(ntests,1);
dfvAL = zeros(ntests,1);
dfvOL = zeros(ntests,1);

for t = 1:ntests
    fprintf('test %d of %d\n',t,ntests);

    AA = var_specrad(randn(nvars,nvars,morder),rho);
    VV = random_covmat(nvars,corrfac);
    X = var_to_tsdata(AA,VV,nobs,ntrials);

    [AO,VO] = tsdata_to_var(X,morder,'OLS');
    [AL,VL] = tsdata_to_var(X,morder,'LWR');

    dfaAO(t) = maxabs(AA-AO);
    dfaAL(t) = maxabs(AA-AL);
    dfaOL(t) = maxabs(AL-AO);

    dfvAO(t) = maxabs(VV-VO);
    dfvAL(t) = maxabs(VV-VL);
    dfvOL(t) = maxabs(VL-VO);
end

fprintf('\nA diff actual - OLS = %.8f   +-   %.8f \n',mean(dfaAO),std(dfaAO));
fprintf('A diff actual - LWR = %.8f   +-   %.8f \n',mean(dfaAL),std(dfaAL));
fprintf('A diff OLS    - LWR = %.8f   +-   %.8f \n',mean(dfaOL),std(dfaOL));

fprintf('\nV diff actual - OLS = %.8f   +-   %.8f \n',mean(dfvAO),std(dfvAO));
fprintf('V diff actual - LWR = %.8f   +-   %.8f \n',mean(dfvAL),std(dfvAL));
fprintf('V diff OLS    - LWR = %.8f   +-   %.8f \n',mean(dfvOL),std(dfvOL));
