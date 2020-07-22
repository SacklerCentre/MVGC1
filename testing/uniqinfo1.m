n = 1000;

a = rand(n,1);
b = rand(n,1);
c = rand(n,1);

a = sqrt((1+b)/2);
c = a;

LHS = (1-a.*a).*(1-b.*b).*(1-c.*c);
RHS = 1 + 2*a.*b.*c - (a.*a + b.*b + c.*c);

violations = LHS > RHS;

fprintf('number of violations = %d\n',nnz(violations));
[a b c violations];
