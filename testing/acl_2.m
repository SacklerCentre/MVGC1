n = 2;
q = 8;

aclags = 40;

%----------------------------

V = eye(n);

B = cat(3,eye(n),randn(n,n,q));

for k = 0:q
    Gk = zeros(2);
    for m = 0:q-k
        Gk = Gk + B(:,:,k+m+1)*B(:,:,m+1)';
    end
    G(:,:,k+1) = Gk;
end
G = cat(3,G,zeros(n,n,aclags-q));

figure(1); clf;
plot_autocov(G);

ptic('*** autocov_to_var\n');
[A1,V1] = autocov_to_var_test(G);
ptoc('*** autocov_to_var took ');

figure(1); clf;
plot_autocov(V1,[],[],[],false,false);
