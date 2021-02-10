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

AA = [reshape(A,n,pn); eye(pn1) zeros(pn1,n)];
H  = var2trfun(A,fres);
PSIGL = zeros(n1,n1,n);
for j = 1:n
    jo = [1:j-1 j+1:n]; % omit j
    PSIGL(:,:,j) = chol(parcov(SIG,jo,j),'lower'); % pre-compute the partial covariances for efficiency
end

for j = 1:n
    jo = [1:j-1 j+1:n]; % omit j

	% Solve the shrunken DARE

	[KJ,SIGj,rep] = var2riss(A,SIG,j,jo);
    if sserror(rep,j), continue; end % check DARE report, bail out on error

	% Calculate reduced SS parameters from shrunken DARE (note: SIGj is the same)

	Cj = reshape(A(jo,:,:),n1,pn);
	Kj = zeros(pn,n1);
	Kj(jo,:) = eye(n1);
	qn = 0;
	for q = 1:p
		Kj(qn+j,:) = KJ(q,:);
		qn = qn+n;
	end
    Bj = ss2itrfun(AA,Cj,Kj,fres);

	% Calculate spectral GC

    for ii = 1:n1
        i  = jo(ii);
        ijo = [1:i-1 i+1:n]; % omit i
        Sj  = SIGj(ii,ii);   % reduced model spectrum is flat!
        LSj = log(Sj);
        for k = 1:h
            Hjk = Bj(ii,:,k)*H(jo,ijo,k)*PSIGL(:,:,i);
            f(i,j,k) = LSj - log(Sj-Hjk*Hjk');
        end
    end
end
