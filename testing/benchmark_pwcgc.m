S = 100;

n = 10;
m = 1000;
p = 20;

regm = 'OLS';

minr = 0.5;
maxr = 0.9;
I    = 5;

%-----------------------------------------

r = linspace(minr,maxr,I)';

t = zeros(I,3);
for i = 1:I;
    fprintf('iter %d of %d : r = %f\n',i,I,r(i));
    
    for s = 1:S
        A = var_specrad(randn(n,n,p),r(i));
        V = randn(n); V = V*V';
        X = var_to_tsdata(A,V,m);
%{        
        tic;
        [AA,VV] = tsdata_to_var(X,p,regm);
        t(i,1) = t(i,1) + toc;
        
        tic;
        G = var_to_autocov(AA,VV);
        t(i,2) = t(i,2) + toc;
        
        tic;
        F = autocov_to_mvgc(G,1,2);
        t(i,3) = t(i,3) + toc;
%}        
        cput = cputime;
        [AA,VV] = tsdata_to_var(X,p,regm);
        t(i,1) = t(i,1) + (cputime-cput);
        
        cput = cputime;
        G = var_to_autocov(AA,VV);
        t(i,2) = t(i,2) + (cputime-cput);
        
        cput = cputime;
        F = autocov_to_mvgc(G,1,2);
        t(i,3) = t(i,3) + (cputime-cput);
        
    end
end
t = 1000*t./S;

plot(r,t);
title(sprintf('n = %d, m = %d, p = %d',n,m,p));
xlabel('spectral radius');
ylabel('CPU time (ms)');
legend('tsdata\_to\_var','var\_to\_autocov','autocov\_to\_gc','Location','NorthWest');
