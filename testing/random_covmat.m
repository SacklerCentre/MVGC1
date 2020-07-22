function V = random_covmat(n,grho)

I = eye(n);
if grho < eps, V = I; return; end
C = randn(n);
S = cov2corr(C*C')-I;
c = grho*sqrt(2/trace(S*S));
V = I + c*S;

% now sqrt(1-det(V)) ~ grho, at least for smallish grho
