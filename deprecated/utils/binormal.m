function X = binormal(m,n,parms)

% X = binormal(m,n,parms)
%
% Return random variates from a superposition of two Gaussian distributions,
% with parms.p the weight on the first distribution.

P = (rand(m,n) < parms.p);
X = P.*(parms.amu+parms.asig*randn(m,n))+(~P).*(parms.bmu+parms.bsig*randn(m,n));

