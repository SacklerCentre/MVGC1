function [x0,i0] = eroot(y,x)

% [x0,i0] = eroot(y,x)
%
% Find values of x0 in x for which y = 0 (linear interpolates). The (lower)
% index corresponding to the zero crossing is returned in i0.

n = length(y);
assert(isvector(y) && isvector(x) && length(x) == n,'y,x must be vectors of the same length');
yp  = y>0;
yp1 = [yp(2:end);yp(end)];
yn  = y<=0;
yn1 = [yn(2:end);yn(end)];
i0  = find((yp&yn1)|(yp1&yn));
if isempty(i0)
    x0 = [];
else
    assert(max(i0) < n,'uh oh');
    x0  = (x(i0).*y(i0+1)-x(i0+1).*y(i0))./(y(i0+1)-y(i0));
end
