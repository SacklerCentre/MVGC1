n = 2;
m = 1000;
mtrunc = 1000;
p = 4;
w = 0.005;
a1 = 20;
a2 = -13;

rho = 0.95;

%-------------------------------------------------

m = m+mtrunc;

rng_seed(s);

AA = randn(n,n,p);
A = zeros(n,n,p,m);
for t = 1:m
    AA(1,1,1) = a1*sin(w*t);
    AA(1,2,1) = a1*sin(w*t);
    AA(2,1,1) = a2*sin(w*(t-100));
    A(:,:,:,t) = var_specrad(AA,rho);
    
end

X = randn(n,m);
for t = p+1:m
    for k = 1:p
        X(:,t) = X(:,t) + A(:,:,k,t)*X(:,t-k);
    end
end

m = m-mtrunc;
X = X(:,end-m+1:end);
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
