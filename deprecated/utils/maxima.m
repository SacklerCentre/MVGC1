function [xmax,imax] = maxima(x,ignorenans)

x = x(:)';
n = length(x);
nans = isnan(x);
if nargin > 1 && ignorenans
    y = x;
    y(nans) = [];
    lmax = false(1,n);
    lmax(~nans) = diff(sign(diff([-Inf y -Inf]))) < 0;
    imax = 1:n;
    imax = imax(lmax);
else
    x(nans) = -Inf;
    imax = 1:n;
    imax = imax(diff(sign(diff([-Inf x -Inf]))) < 0);
end
xmax = x(imax);
