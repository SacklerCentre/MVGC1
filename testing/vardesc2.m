
N         = 100;
n         = 9;
n1        = 4;
p         = 3;
q         = 10;
rho       = 0.9;    % spectral radius
grho      = 0.5;    % residuals correlation factor

aseed     = 193873;
vseed     = 912982;
tseed     = 0;

%-------------------------------------------------------------------------------

s = rng_seed(aseed);
A = var_specrad(randn(n,n,p),rho);
rng_restore(s);

s = rng_seed(vseed);
V = random_covmat(n,grho);
rng_restore(s);

s = rng_seed(tseed);
[FF,TT,res] = vardesc(A,V,n1,N);
rng_restore(s);

fprintf('res = '); disp(res);

fprintf('FF = %g\n',FF);
