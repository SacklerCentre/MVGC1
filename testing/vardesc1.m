regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)

acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)

tstat     = '';     % statistical test for MVGC:  'chi2' for Geweke's chi2 test (default) or'F' for Granger's F-test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

m         = 100000;  % number of observations

n         = 5;
n1        = 2;
p         = 3;
q         = 10;
rho       = 0.9;    % spectral radius
grho      = 0.5;    % residuals correlation factor

aseed     = 193873;
vseed     = 912982;
tseed     = 0;
xseed     = 0;

normt = 1;
vfcheck = 0;
moest = 0;

%-------------------------------------------------------------------------------

seed      = 0;      % random seed (0 for unseeded)

rng_seed(seed);

s = rng_seed(aseed);
A = var_specrad(randn(n,n,p),rho);
rng_restore(s);

s = rng_seed(vseed);
V = random_covmat(n,grho);
rng_restore(s);

s = rng_seed(tseed);
T = randn(n1,n);
rng_restore(s);

if normt
    T = chol(inv(T*V*T'))*T; % so T*V*T' = I
%   fprintf('logdet(TVT'') = %g\n',log(det(T*V*T')));
end

s = rng_seed(xseed);
X = var_to_tsdata(A,V,m);
rng_restore(s);

TX = T*X;

if moest
    momax = 20;
    
    [AIC,BIC] = tsdata_to_infocrit(TX,momax);
    [~,bmo_AIC] = min(AIC);
    [~,bmo_BIC] = min(BIC);

    figure(1); clf;
    plot((1:momax)',[AIC BIC]);
    legend('AIC','BIC');

    fprintf('\nbest model order (AIC) = %d\n',bmo_AIC);
    fprintf('best model order (BIC) = %d\n',bmo_BIC);
    return
end

VF = T*V*T';

if vfcheck
    [~,VF1] = lregress(TX,X,q);
    fprintf('VF =\n'); disp(VF);
    fprintf('VF1 =\n'); disp(VF1);
end

[~,VR] = tsdata_to_var(TX,q,regmode);
assert(~isbad(VR),'reduced VAR estimation failed');

F = log(det(VR))-log(det(VF));

fprintf('F = %g\n',F);
