
regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

acmaxlags = 2000;   % maximum autocovariance lags (empty for automatic calculation)

n = 10;
p = 5;

x = 1;

rho = 0.99;

A = var_specrad(1-2*rand(n,n,p),rho);

% Residuals covariance matrix.

SIG = eye(n);

ptic('*** var_to_autocov... ');
[G,info] = var_to_autocov(A,SIG,acmaxlags);
ptoc;

var_info(info,true); % report results (and bail out on error)

GR = G(x,x,:);

ptic('*** var_to_autocov... ');
[AR,SIGR] = autocov_to_var(GR);
ptoc;

maxa = max(abs(A(:)));

figure(1); clf;
plot_varcoeffs(AR(:,:,1:40)/maxa);
ylim([-1 1]);
