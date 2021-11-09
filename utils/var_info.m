function info = var_info(A,V,report)

if nargin < 3 || isempty(report), report = 1; end % default: print out report

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
[n1,n2]  = size(V);
assert(n1 == n && n2 == n,'Residuals covariance matrix must be square and match VAR coefficients matrix');
pn1 = (p-1)*n;

info.error = uint32(0);

info.observ = n;
info.morder = p;

% calculate spectral radius

info.rho = max(abs(eig([reshape(A,n,p*n); eye(pn1) zeros(pn1,n)],'nobalance'))); % v2.0 - don't balance!

info.acdec = ceil(0.5*log(eps)/log(info.rho)); % so that autocov decays to < sqrt(eps), (probably ~ 1.5e-8)

if maxabs(triu(V,1)-triu(V',1)) > eps
    info.sigspd = 1; % not symmetric
else
    [~,cholp] = chol(V,'lower');
    if cholp > 0
        info.sigspd = 2; % symmetric, but not positive definite
    else
        info.sigspd = 0; % symmetric, positive definite
    end
end
info.mii = multiinfo(V);        % multi-information (generalised correlation)
info.umii = multiinfo(n,true);  % multi-information for uniform random n x n correlation matrix, for comparison

%{
if rand < 0.1 % test intermittent errors
    info.mii = -0.5;
end
%}

rhotol = sqrt(eps);

if     info.rho > 1+rhotol, disp('ex'); info.error = bitset(info.error,1); % explosive
elseif info.rho > 1-rhotol, disp('ur'); info.error = bitset(info.error,2); % unit root
end

if     info.sigspd == 1,     info.error = bitset(info.error,5); % not symmetric
elseif info.sigspd == 2,     info.error = bitset(info.error,6); % not positive definite
end

if     info.mii < 0,         info.error = bitset(info.error,7); % negative
end

if report == 1 % print out report

    fprintf('\nVAR info:\n');

    fprintf('    variables         = %d\n',info.observ);

    fprintf('    model order       = %d\n',info.morder);

    fprintf('    AR spectral norm  = %.4f',info.rho);
    if      bitget(info.error,1), fprintf(2,'    ERROR: unstable (explosive)\n');
    elseif  bitget(info.error,2), fprintf(2,'    ERROR: unstable (unit root)\n');
    else    fprintf('    stable (autocorrelation decay ~ %d)\n',info.acdec);
    end

    fprintf('    residuals covariance matrix');
    if     bitget(info.error,5), fprintf(2,'     ERROR: not symmetric\n');
    elseif bitget(info.error,6), fprintf(2,'     ERROR: not positive definite\n');
    elseif bitget(info.error,7), fprintf(2,'     ERROR: multi-information negative\n');
    else   fprintf('   symmetric, pos. def. (mii = %.4f, uniform = %.4f)\n',info.mii,info.umii);
    end

    fprintf('\n');

elseif report > 1 % format error message(s) string

    if ~info.error, info.errmsg = ''; return; end % no errors to report

    info.nerrors = nnz(bitget(info.error,1:8)); % number of errors

    if info.nerrors > 1
        info.errmsg = 'VAR ERRORS';
    else
        info.errmsg = 'VAR ERROR';
    end

    if      bitget(info.error,1), info.errmsg = [info.errmsg sprintf(': AR spectral norm = %.6f - unstable (explosive)',info.rho)];
    elseif  bitget(info.error,2), info.errmsg = [info.errmsg sprintf(': AR spectral norm = %.6f - unstable (unit root)',info.rho)];
    end

    if     bitget(info.error,5), info.errmsg = [info.errmsg ': res. cov. matrix not symmetric'];
    elseif bitget(info.error,6), info.errmsg = [info.errmsg ': res. cov. matrix not positive definite'];
    end

    if     bitget(info.error,7), info.errmsg = [info.errmsg sprintf(': multi-information = %.6f - negative',info.mii)];
    end

end
