function [C,p] = normalise_chol(C,cfac,sdev)

% C = normalise_chol(C,cfac,sdev)
%
% Normalise a Cholesky matrix C so that the standard deviations of the
% corresponding covariance matrix C*C' are all equal to sdev. Also apply
% the factor cfac (should be <= 1) to correlations.

n = size(C,1);

SIG = C*C';
T = sdev*diag(1./sqrt(diag(SIG)));
SIG = T*SIG*T;
SIG = cfac*SIG;
SIG(1:n+1:n*n) = ones(n,1);
[C,p] = chol(SIG,'lower');

