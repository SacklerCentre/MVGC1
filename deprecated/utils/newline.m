function newline(n)

if nargin < 1
    n = 1;
end

fprintf(repmat('\n',1,n));
