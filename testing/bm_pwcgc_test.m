S = 10;

n = 10;
m = 1000;
N = 1;
p = 3;

regm = 'OLS';

r = 0.9;

%-----------------------------------------

pp = 1:10;
I = length(pp);

t = zeros(2,I);
for i = 1:I;
    fprintf('iter %d of %d\n',i,I);
    p = pp(i);
    [t(1,i),t(2,i)] = bm_pwcgc(S,n,m,N,p,r,regm);
end

plot(pp,t);
legend('MVGC','GCCA');
