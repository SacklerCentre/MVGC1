S = 10000;

n = 2;
m = 1000;
N = 1;
p = 1;

regm = 'OLS';
acmaxl = 1000;

a    = 0.8;
b    = 0.9;

cmax = 4;
I    = 21;

%rng_seed(28272);

%-------------------------------------------------------------------------------

c = linspace(0,cmax,I)';

A = [a 0; 0 b];
V = eye(n);

F1c = zeros(S,I);
F1n = zeros(S,I);
F2c = zeros(S,I);
F2n = zeros(S,I);

t = zeros(1,I);
for i = 1:I;
    fprintf('iteration %d of %d : c = %f ',i,I,c(i));
    
    A(1,2) = c(i);
    
    for s = 1:S
        X1       = var_to_tsdata(A,V,m);
        [A1,V1]  = tsdata_to_var(X1,p,regm);
        [G,info] = var_to_autocov(A1,V1,acmaxl);
        while info.error
            fprintf('*');
            X1       = var_to_tsdata(A,V,m);
            [A1,V1]  = tsdata_to_var(X1,p,regm);
            [G,info] = var_to_autocov(A1,V1,acmaxl);
        end
        FF       = autocov_to_pwcgc(G);
        F1c(s,i) = FF(1,2);
        F1n(s,i) = FF(2,1);

        X2        = var_to_tsdata(A,V,m);
        [FF,~,V2] = GCCA_tsdata_to_pwcgc(X2,p,regm);
        while isbad(V2)
            fprintf('#');
            X2        = var_to_tsdata(A,V,m);
            [FF,~,V2] = GCCA_tsdata_to_pwcgc(X2,p,regm);
        end
        F2c(s,i) = FF(1,2);
        F2n(s,i) = FF(2,1);
    end
    fprintf('\n');
end

D = 1+b^2+c.^2;
Finf = log((D + sqrt(D.^2 - 4*b^2))/2);
Finf(Finf < 0) = 0;

b0 = 1/(1-b^2);
c0 = (b*c./(1-a*b))*b0;
a0 = (1+(c.^2)*b0+2*a*c.*c0)/(1-a^2);
a1 = a*a0+c.*c0;
Fone = log((a0.^2-a1.^2)./a0);
Fone(Fone < 0) = 0;

save(sprintf('accu_test_c__m_%d_a_%.2f_b_%.2f_cmax%.2f_I_%d_S_%d.mat',m,a,b,cmax,I,S));
