S = 100;

n = 10;
m = 1000;
N = 1;
p = 5;

regm = 'LWR';

r = 0.9;

%-----------------------------------------

pp = 1:20;
I = length(pp);

t = zeros(1,I);
for i = 1:I;
    fprintf('iter %d of %d\n',i,I);
    p = pp(i);
    t(i) = bm_tsdata_to_var(S,n,m,N,p,r,regm);
end

plot(pp,t);
