function [t1,t2] = bm_pwcgc(S,n,m,N,p,r,regm)

X = zeros(n,m,N,S);
for s = 1:S
    A = var_specrad(randn(n,n,p),r);
    V = randn(n); V = V*V';
    X(:,:,:,s) = var_to_tsdata(A,V,m,N);
end

tic;
for s = 1:S
    [A1,V1,E1] = tsdata_to_var(X(:,:,:,s),p,regm);
    G = var_to_autocov(A1,V1);
    F1 = autocov_to_pwcgc(G);
end
t1 = toc/S;

tic;
for s = 1:S
    [F2,A2,V2,E2] = GCCA_tsdata_to_pwcgc(X(:,:,:,s),p,regm);
end
t2 = toc/S;
