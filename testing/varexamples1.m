n = 2;
m = 2000;
p = 4;

w1 = 0.05;
w2 = 0.010;
a1 = 2;
a2 = -1;

dt1 = 0.01;
dt2 = -0.002;

rho = 0.9;

%-------------------------------------------------

rng_seed(s);

A = var_specrad(randn(n,n,p),rho);

X = var_to_tsdata(A,eye(n),m);

T = 1:m;

X(1,:) = X(1,:) + dt1*T + a1*sin(w1*T);
X(2,:) = X(2,:) + dt2*T + a1*sin(w2*(T-100));

figure(1); clf; plot(T,X);

momax = 70;
[AIC,BIC] = tsdata_to_infocrit(X,momax,'LWR');
[~,bmo_AIC] = min(AIC);
[~,bmo_BIC] = min(BIC);

fprintf('\nbest model order (AIC) = %d\n',bmo_AIC);
fprintf('best model order (BIC) = %d\n',bmo_BIC);

figure(2); clf;
plot_tsdata([AIC BIC]',{'AIC','BIC'});


[AAA,VVV] = tsdata_to_var(X,bmo_AIC);
[G,info] = var_to_autocov(AAA,VVV);
var_info(info);
