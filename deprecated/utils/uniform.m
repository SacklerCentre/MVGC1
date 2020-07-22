function X = uniform(m,n,parms)

% X = uniform(m,n,parms)
%
% Return random variates fro a uniform distribution.

X = parms.lo+(parms.hi-parms.lo)*rand(m,n);

