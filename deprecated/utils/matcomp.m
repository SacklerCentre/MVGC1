function s = matcomp(X,Y,longfmt)

assert(isequal(size(X),size(Y)),'matrices don''t match');

if nargin < 3 || isempty(longfmt)
    longfmt = false;
end

diff = abs(X(:)-Y(:));

dmax  = max(diff);
dmean = mean(diff);
dstd  = std(diff);

if longfmt
    s = sprintf('mean = %.16e %c %.16e, max = %.16e, ',dmean,char(177),dstd,dmax);
else
    s = sprintf('mean = %e %c %e, max = %e',dmean,char(177),dstd,dmax);
end
