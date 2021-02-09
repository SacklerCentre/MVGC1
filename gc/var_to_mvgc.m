function F = var_to_mvgc(A,SIG,x,y)

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(SIG);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match VAR coefficients matrix');

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(length(unique([x y])) == length([x y]),'x and y indices must be unique and non-overlapping');
assert(all(x >=1 & x <= n),'Some x indices out of range');
assert(all(y >=1 & y <= n),'Some y indices out of range');

z  = 1:n; z([x y]) = []; % indices of other variables (to condition out)
xz = [x z];               % indices of variables in reduced model (omit source variables)

nx = length(x);
ny = length(y);
xr = 1:nx; % indices of x in reduced model

F = NaN;
[~,SIGR,rep] = var2riss(A,SIG,y,xz); % residuals covariance matrix of reduced model
if sserror(rep), return; end         % check DARE report, bail out on error
F = logdet(SIGR(xr,xr)) - logdet(SIG(x,x));
