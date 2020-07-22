function [X,E] = arnonlin(a,b,c,e,E,trunc,u,v)

m = size(E,2);

if nargin >= 8  && ~isempty(u) && ~isempty(v)
    E = [u 0; 0 v]*E; 
end

X = E;
for t = 2:m
    X(1,t) = X(1,t) + a*X(1,t-1) + c*real(X(2,t-1)^e);
    X(2,t) = X(2,t) + b*X(2,t-1);
end

if nargin >= 6 && ~isempty(trunc) && trunc > 0
    X = X(:,trunc+1:m);
    if nargout > 1
        E = E(:,trunc+1:m);
    end
end
