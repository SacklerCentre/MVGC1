function f = var_to_spwcgc(A,SIG,fres)

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2] = size(SIG);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square, and match coefficients matrix');

n1  = n-1;
pn  = p*n;
pn1 = pn-n;

h = fres+1;
f = nan(n,n,h);

% for efficiency we pre-compute some stuff

AR = [reshape(A,n,pn); eye(pn1) zeros(pn1,n)];
H  = var2trfun(A,fres);
PSIGL = zeros(n1,n1,n);
for y = 1:n
    r = [1:y-1 y+1:n]; % omit y
    PSIGL(:,:,y) = chol(parcov(SIG,r,y),'lower'); % pre-compute the partial covariances for efficiency
end

for y = 1:n
    r = [1:y-1 y+1:n]; % omit y

	% Solve the shrunken DARE

	[KT,VR,rep] = var2riss(A,SIG,y,r);
    if sserror(rep,y), continue; end % check DARE report, bail out on error

	% Calculate reduced SS parameters from shrunken DARE (note: VR is the same)

	CR = reshape(A(r,:,:),n1,pn);
	KR = zeros(pn,n1);
	KR(r,:) = eye(n1);
	kn = 0;
	for k1 = 1:p
		KR(kn+y,:) = KT(k1,:);
		kn = kn+n;
	end
    BR = ss2itrfun(AR,CR,KR,fres);

	% Calculate spectral GC

    for xr = 1:n1
        x  = r(xr);
        w = [1:x-1 x+1:n];  % omit x

        SR  = VR(xr,xr); % reduced model spectrum is flat!
        LSR = log(SR);

        for k = 1:h
            HR = BR(xr,:,k)*H(r,w,k)*PSIGL(:,:,x);
            f(x,y,k) = LSR - log(SR-HR*HR');
        end
    end
end
