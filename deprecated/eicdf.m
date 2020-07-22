function x = eicdf(y,p)

% x = eicdf(y,p)
%
% Return value of inverse empirical CDF for sample distributions in columns of y,
% evaluated at probabilities p (which must be scalar, or a vector of length
% equal to the number of columns of y).

cols = size(y,2);
if cols > 1
    x = zeros(1,cols);
    if isscalar(p)
        for j=1:cols
            x(j) = eicdf(y(:,j),p);
        end
    else
        assert(isvector(p) && length(p) == cols,'p doesn''t match y');
        for j=1:cols
            x(j) = eicdf(y(:,j),p(j));
        end
    end
    return
end

y = y(~isnan(y)); % remove NaNs

[f,x] = ecdf(y);
x = eroot(f-p,x);
if     length(x) < 1
    fprintf(2,'WARNING: no inverse found\n');
elseif length(x) > 1
    fprintf(2,'WARNING: multiple inverses found!?\n');
end
