function [xmin,imin] = minima(x,ignorenans)

x = x(:)';
n = length(x);
nans = isnan(x);
if nargin > 1 && ignorenans
    y = x;
    y(nans) = [];
    lmin = false(1,n);
    lmin(~nans) = diff(sign(diff([Inf y Inf]))) > 0;
    imin = 1:n;
    imin = imin(lmin);
else
    x(nans) = Inf;
    imin = 1:n;
    imin = imin(diff(sign(diff([Inf x Inf]))) > 0);
end
xmin = x(imin);
