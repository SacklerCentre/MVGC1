seed = 1239123;
%seed = 0;

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)

n = 3;
p = 8;

x = 1;

rho = 0.99;

%-------------------------------------------------------------------------------

rng_seed(seed);


A = var_specrad(1-2*rand(n,n,p),rho);

% Residuals covariance matrix.

SIG = eye(n);

ptic('*** var_to_autocov... ');
[G,info] = var_to_autocov(A,SIG,acmaxlags);
ptoc;

var_info(info,true); % report results (and bail out on error)

GR = G(x,x,:);

ptic('*** var_to_autocov... ');
[AR,SIGR] = autocov_to_var_diag(GR);
ptoc;

q = size(AR,3);

fprintf('\neffective mo = %d\n',q);

maxa = max(abs(AR(:)));

figure(1); clf;
plot_varcoeffs(AR/maxa);
ylim([-1 1]);

rdiff = nan(1,q);
rdmin = nan(1,q);
rdpeak = false(1,q);
for k = 3:q
    %rdiff(k) = abs(norm(SIGR(:,:,k)-SIGR(:,:,k-1)))/norm(SIGR(:,:,k-1));
    rdiff(k) = abs(det(SIGR(:,:,k))-det(SIGR(:,:,k-1)))/det(SIGR(:,:,k-1));
    rdmin(k) = rdiff(k)-rdiff(k-1);
    rdpeak(k-1) = rdmin(k) < 0 &&  rdmin(k-1) > 0;
end

mo = 1:k;
mop = mo(rdpeak);
rdiffp = rdiff(rdpeak);

figure(2); clf;
semilogy(mo,rdiff,mop,rdiffp);
ylim([1e-20 1]);
