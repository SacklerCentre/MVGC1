function X = biuniform(m,n,parms)

% X = biuniform(m,n,parms)
%
% Return random variates fro a superposition of two uniform distributions,
% with parms.p the weight on the first distribution.

P = (rand(m,n) < parms.p);
X = P.*(parms.alo+(parms.ahi-parms.alo)*rand(m,n))+(~P).*(parms.blo+(parms.bhi-parms.blo)*rand(m,n));

