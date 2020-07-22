function ptime(s1,s2)

if nargin <2 || isempty(s2)
    fprintf(1,[s1 datestr(now)]);
else
    fprintf(1,[s1 datestr(now) s2]);
end
