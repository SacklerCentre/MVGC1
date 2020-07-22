function t = bm_tsdata_to_var(S,n,m,N,p,r,regm)

for s = 1:S
    A = var_specrad(randn(n,n,p),r);
    V = randn(n); V = V*V';
    X(:,:,:,s) = var_to_tsdata(A,V,m,N);
end

tic;
for s = 1:S
    [AA,VV,EE] = tsdata_to_var(X(:,:,:,s),p,regm);
end
t = toc/S;
