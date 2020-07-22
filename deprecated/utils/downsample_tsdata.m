function XD = downsample_tsdata(X,kdt,ds)

% XD = downsample_tsdata(X,k)
%
% XD = downsample_tsdata(X,dt,ds)
%
% Downsample time series data X.
%
% For the first form the (integer) downsample factor k is supplied.
%
% For the second form, the downsample factor is calculated from the step sizes
% dt of the original and ds of the downsampled time series.

if nargin < 3
    assert(kdt == floor(kdt),'downsample factor must be an integer');
    k = kdt;
else
    dt = kdt;         % dt = time step (raw)
    k = round(ds/dt); % ds = time step (downsampled)
end
assert(k > 0,         'downsample factor too small');
assert(k <= size(X,2),'downsample factor too large');

XD = X(:,1:k:end);
