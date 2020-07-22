function X = normal(m,n,parms)

%  X = normal(m,n,parms)
%
% Return random variates from a normal distributions.

X =  parms.mu+parms.sig*randn(m,n);

