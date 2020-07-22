function [xm,xlo,xhi] = mconfint(x,alpha)

% [xm,xlo,xhi] = mconfint(x)
%
% Return mean and confidence intervals for data in x

xm  = mean(x);
xlo = eicdf(x,alpha);
xhi = eicdf(x,1-alpha);
