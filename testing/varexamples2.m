n = 2;
m = 2000;
p = 3;

rho = 1;

%-------------------------------------------------

rng_seed(s);

A = var_specrad(randn(n,n,p),rho);

X = var_to_tsdata(A,eye(n),m,1,100);

T = 1:m;

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
