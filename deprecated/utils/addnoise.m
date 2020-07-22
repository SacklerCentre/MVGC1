function Y = addnoise(X,noise,scaled)

% Y = add_snoise(X,noise)
%
% Add Gaussian white noise to X, optionally scaled by std dev.

if nargin < 3 || isempty(noise), scaled = false; end % default: no scaling

[n,m,N] = size(X);

if isscalar(noise)
    noise = noise*ones(1,n);
else
    assert(isvector(noise) && length(noise) == n,'noise must be a scalar or a vector of same length as number of variables');
end
assert(all(noise >= 0) ,'noise must be non-negative');
noise = noise(:); % make col vector

M = N*m;
X = reshape(X,n,M);

if scaled
    Y = X + ((noise.*std(X,[],2))*ones(1,M)).*randn(n,M);
else
    Y = X + (noise*ones(1,M)).*randn(n,M);
end

Y = reshape(Y,n,m,N);
