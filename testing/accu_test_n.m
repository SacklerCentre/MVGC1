S = 100;

n = 2;
m = 1000;
N = 1;

regm = 'OLS';
acmaxl = 1000;

a    = 0.8;
b    = 0.9;
c    = 1;

numax = 4;
I     = 21;

icregm = 'LWR';
momax = 40;

%rng_seed(28272);

%-------------------------------------------------------------------------------

nu = linspace(0,numax,I)';

A = [a c; 0 b];
V = eye(n);

mo1 = zeros(S,I);
F1c = zeros(S,I);
F1n = zeros(S,I);

mo2 = zeros(S,I);
F2c = zeros(S,I);
F2n = zeros(S,I);

t = zeros(1,I);
for i = 1:I;
    fprintf('iteration %d of %d : nu = %f',i,I,nu(i));
    
    for s = 1:S
        X1       = var_to_tsdata(A,V,m) + nu(i)*randn(n,m);
        [~,BIC]  = tsdata_to_infocrit(X1,momax,icregm,false);
        [~,pBIC] = min(BIC);
        mo1(s,i) = pBIC;
        [A1,V1]  = tsdata_to_var(X1,mo1(s,i),regm);
        [G,info] = var_to_autocov(A1,V1,acmaxl);
        while info.error
            fprintf('*');
            X1       = var_to_tsdata(A,V,m) + nu(i)*randn(n,m);
            [~,BIC]  = tsdata_to_infocrit(X1,momax,icregm,false);
            [~,pBIC] = min(BIC);
            mo1(s,i) = pBIC;
            [A1,V1]  = tsdata_to_var(X1,mo1(s,i),regm);
            [G,info] = var_to_autocov(A1,V1,acmaxl);
        end
        FF       = autocov_to_pwcgc(G);
        F1c(s,i) = FF(1,2);
        F1n(s,i) = FF(2,1);

        X2        = var_to_tsdata(A,V,m) + nu(i)*randn(n,m);
        [~,BIC]   = tsdata_to_infocrit(X2,momax,icregm,false);
        [~,pBIC]  = min(BIC);
        mo2(s,i)  = pBIC;
        [FF,~,V2] = GCCA_tsdata_to_pwcgc(X2,mo2(s,i),regm);
        while isbad(V2)
            fprintf('#');
            X2        = var_to_tsdata(A,V,m) + nu(i)*randn(n,m);
            [~,BIC]   = tsdata_to_infocrit(X2,momax,icregm,false);
            [~,pBIC]  = min(BIC);
            mo2(s,i)  = pBIC;
            [FF,~,V2] = GCCA_tsdata_to_pwcgc(X2,mo2(s,i),regm);
        end
        F2c(s,i) = FF(1,2);
        F2n(s,i) = FF(2,1);
    end
    fprintf('\n');
end

D = 1+b^2+c.^2;
Finf = log((D + sqrt(D.^2 - 4*b^2))/2);
Finf = repmat(Finf,I,1);

b0 = 1/(1-b^2);
c0 = (b*c./(1-a*b))*b0;
a0 = (1+(c.^2)*b0+2*a*c.*c0)/(1-a^2);
a1 = a*a0+c.*c0;
Fone = log((a0.^2-a1.^2)./a0);
Fone = repmat(Fone,I,1);

save(sprintf('accu_test_n__m_%d_a_%.2f_b_%.2f_c%.2f_numax_%.2f_I_%d_S_%d.mat',m,a,b,c,numax,I,S));
