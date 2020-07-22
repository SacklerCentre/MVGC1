function [X,E] = genvma(B,E,trunc)

% [X,E] = genvma(B,E,trunc)
%
% Generate vector moving average process (VMA) with coefficients B and residuals
% E - implements:
%
%     X(t) = E(t) + sum(k = 1:p) B(k)*E(t-k)
%
% Initial values for the outputs X(t) are set to the residuals E(t). Optionally
% truncate the first 'trunc' values.

if nargin < 3 || isempty(trunc); trunc = 0; end

[n,m] = size(E);
[n1,n2,p] = size(B);
assert(n1 == n2, 'coefficients matrix not square');
assert(n1 == n,  'coefficients matrix doesn''t match residuals');
assert(trunc >= 0 && trunc < m,'bad truncation');

m1 = m-p;
p1 = p+1;

% stack 0 to p lags of residuals

EE = zeros(n,p1,m1);
for k = 0:p
    EE(:,k+1,:) = E(:,p1-k:m-k); % concatenate k-lagged residuals
end
EE = reshape(EE,n*p1,m1); % stack lagged residuals

X = [E(:,1:p) [eye(n) B(:,:)]*EE]; % the VMA

if trunc > 0
    X = X(:,trunc+1:m);
    if nargout > 1
        E = E(:,trunc+1:m);
    end
end
