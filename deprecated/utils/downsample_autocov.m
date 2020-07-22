function GD = downsample_autocov(G,q,kdt,ds)

% GD = downsample_autocov(G,q,k)
%
% GD = downsample_autocov(G,q,dt,ds)
%
% Downsample autocovariance sequence G, up to q (downsampled) lags
%
% For the first form the (integer) downsample factor k is supplied.
%
% For the second form, the downsample factor is calculated from the step sizes
% dt of the original and ds of the downsampled autocovariance sequence.

if nargin < 4
    assert(kdt == floor(kdt),'downsample factor must be an integer');
    k = kdt;
else
    dt = kdt;         % dt = time step (raw)
    k = round(ds/dt); % ds = time step (downsampled)
end
assert(k > 0,             'downsample factor too small');
assert(1+k*q <= size(G,3),'downsample factor too large');

GD = G(:,:,1+(k*(0:q)));
