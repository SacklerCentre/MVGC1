function [fpos,fneg,nfpos,nfneg] = sigerror(C,S)

% [fpos,fneg,nfpos,nfneg] = sigerror(C,S)
%
% Return matrix of false positives/negatives and counts from actual connection
% matrix C and significance matrix S. Diagonals (self-connections) ignored.

n = size(C,1);
assert(size(C,2) == n,'connection matrix not square');
assert(isequal(size(S),size(C)),'matrices don''t match');

d = 1:n+1:n*n;

C(d) = 0; C = logical(C);
S(d) = 0; S = logical(S);

fpos = (~C & S)+0; % make numeric
fneg = (C & ~S)+0; % make numeric

nfpos = sum(fpos(:));
nfneg = sum(fneg(:));

fpos(d) = NaN; % NaNs on diagonal
fneg(d) = NaN; % NaNs on diagonal

