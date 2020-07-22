function q = aclags(rho,acdectol)

if nargin < 2 || isempty(acdectol), acdectol  = 1e-8; end % autocovariance decay tolerance

q = ceil(log(acdectol)/log(rho)); % minimum lags to achieve specified tolerance
