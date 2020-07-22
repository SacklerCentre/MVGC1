function [X,n] = xoutl(X,sdminfac,sdmaxfac)

% [X,n] = xoutl(X,sdminfac,sdmaxfac)
%
% Eliminate outliers of vector X, defined as less than mean minus sdminfac
% or greater than mean plus sdmaxfac times the standard deviation of X. If
% sdmaxfac not supplied it is taken as equal to sdminfac. Also return new
% length of vector X.

assert(isvector(X));

xmean = mean(X);
xsdev = std(X);

xmin = xmean - sdminfac*xsdev;
if nargin == 3
    xmax = xmean + sdmaxfac*xsdev;
else
    xmax = xmean + sdminfac*xsdev;
end

X(X < xmin | X > xmax) = [];

n = length(X);
