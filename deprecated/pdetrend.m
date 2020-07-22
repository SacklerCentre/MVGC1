function [Y,P,p] = pdetrend(X,pdeg,x)

% [X,P] = pdetrend(X,pdeg,x)
%
% Polynomial detrend of degree pdeg (adapted from Matlab `polyfit' and
% `polyval'). Default is linear detrend (NOTE: demeans as well).

if isempty(pdeg) % do nothing
    Y = X;
    P = [];
    return
end

[n,m,N] = size(X);

if nargin < 2 || isempty(pdeg)
    pdeg = 1; % linear detrend (default) (note that pdeg = 0 simply demeans)
end

if nargin < 3 || isempty(x)
    x = 1:m;
end

mu  = mean(x,2);
sig = std(x,[],2);
x = (x - mu)/sig; % normalise

d = pdeg+1;
V = zeros(d,m);  % Vandermonde matrix
V(d,:) = ones(1,m);
for j = pdeg:-1:1
    V(j,:) = x.*V(j+1,:);
end
p = mean(X,3)/V; % mean over trials

P = zeros(n,m);
for i = 1:n
    P(i,:) = p(i,1)*ones(1,m);
    for j=2:d
        P(i,:) = x.*P(i,:) + p(i,j);
    end
end

Y = zeros(size(X));
for r = 1:N
    Y(:,:,r) = X(:,:,r)-P;
end
