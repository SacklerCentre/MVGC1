function [FF,TT,res] = vardesc(A,SIG,n1,nsamps,acmaxlags,acdectol,dlyap_alg,maxiters,maxrelerr)

if nargin < 5, acmaxlags = []; end
if nargin < 6, acdectol  = []; end
if nargin < 7, dlyap_alg = []; end
if nargin < 8, maxiters  = []; end
if nargin < 9, maxrelerr = []; end

[n,nn,~] = size(A);
assert(nn == n,'VAR coefficients matrix has bad shape');

[nn1,nn2] = size(SIG);
assert(nn1 == nn2,'residuals covariance matrix not square');
assert(nn1 == n  ,'residuals covariance matrix doesn''t match VAR coefficients matrix');

[G,res] = var_to_autocov(A,SIG,acmaxlags,acdectol,dlyap_alg,maxiters,maxrelerr);
if res.error, return; end
q = size(G,3);
TG = zeros(n1,n1,q);

FF = realmax;
TT = 0;
for s = 1:nsamps
    
    T = randn(n1,n);
    T = chol(inv(T*SIG*T'))*T; % so T*SIG*T' = I

    for k = 1:q
        TG(:,:,k) = T*G(:,:,k)*T';
    end

    [~,SIGR] = autocov_to_var(TG);

    F = log(det(SIGR));

    fprintf('F = %g\n',F);

    if F < FF
        FF = F;
        TT = T;
    end
end
