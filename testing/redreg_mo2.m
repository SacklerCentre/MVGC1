seed = 1239123;
seed = 0;

acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)

n = 10;
p = 8;

x = 1:2;

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
[AR,SIGR,AA] = autocov_to_var_diag(GR);
ptoc;

q = size(AR,3);

fprintf('\neffective mo = %d\n',q);

rdiff = nan(1,q);
for k = 2:q
    rdiff(k) = abs(norm(SIGR(:,:,k)-SIGR(:,:,k-1)))/norm(SIGR(:,:,k-1));
end

a = nan(1,q);
for k = 1:q
    a(k) = abs(norm(AA(:,:,k)));
end

mo = 1:k;

figure(1); clf;
semilogy(mo,a,mo,rdiff);
