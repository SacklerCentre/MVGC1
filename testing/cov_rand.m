function V = cov_rand(n,icfac,corrmat)

% Generate random positive-definite covariance matrix
%
% n       - observation variable dimension
% icfac   - innovations correlation factor: a positive integer (as icfac gets
%           bigger, innovation correlations get smaller), or set to Inf for zero
%           correlation
% corrmat - flag: return correlation matrix (default)
%
% V       - covariance matrix
%
% V is generated as the sample covariance matrix of an n-dimensional
% uncorrelated Gaussian white noise sequence of length icfac*n. This ensures
% positive-definiteness. As ifac -> Inf, correlation between innovations -> 0.
% Setting icfac = Inf yields a diagonal V (no correlation). If the corrmat flag
% is set, V is normalised to variance 1; i.e. V is returned as a correlation
% matrix (the identity matrix in the case icfac = Inf).

assert(isscalar(icfac) && isnumeric(icfac) && icfac == floor(icfac) && icfac >= 1);

if nargin< 3 || isempty(corrmat), corrmat = true; end

if isinf(icfac) % diagonal covariance matrix
    if corrmat
        V = eye(n);
    else
        SQRTV = randn(n,1);
        V = diag(SQRTV.^2);
    end
else
    SQRTV = randn(n,icfac*n)/sqrt(icfac*n); % sample covariance of n x n white noise of length icfac*n
    V = SQRTV*SQRTV';
    if corrmat
        D = diag(1./sqrt(diag(V)));
        V = D*V*D; % correlation matrix
    end
end

