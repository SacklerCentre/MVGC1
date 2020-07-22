function F = tsdata_to_mvgc_pwc_permtest_fandl(U,p,nperms,qmax,acdecayfac)

% F = tsdata_to_mvgc_pwc_permtest(U,p,nperms,bsize,qmax,acdecayfac)
%
% For each pairwise conditional MVGC from time series data in U, return nperms
% permutation samples in F. p specifies the model order. First index of F
% is permutation number, second is 'to' variable, third is 'from' variable.
%
% If the block size bsize is a non-negative integer, a "naive" block permutation
% method is used, where blocks of specified size of the predictor variable y are
% permuted. If bsize = 0 (default) the recommended block size bsize = p is used,
% since this might be expected to best preserve the (non-causal) correlation
% structure in the data. If bsize is negative, the more correct (but much
% slower) method of Freedman and Lane is used (experimental) - see routine
% 'freedman_lane'.
%
% Remaining parameters are passed to the 'var_to_autocov' routine.

if nargin <  4, qmax       = []; end
if nargin <  5, acdecayfac = []; end

[n,m,N] = size(U);
assert(N == 1,'multi-trial version not implemented');
assert(m > p,'too many lags');

F = nan(nperms,n,n);

for j=1:n;
    jo  = [1:j-1 j+1:n]; % omit j

    Uj = U(:,p+1:m);

    for ii=1:n-1;
        i = jo(ii);

        fprintf('from node %d -> %d ... ',j,i);

        [AX,XZ,EE] = freedman_lane(U,p,i,j);

        nbad = 0;
        for t = 1:nperms

            %fprintf('from node %d -> %d : permutation %d of %d',j,i,t,nperms);

            Uj(i,:) = genvar(AX,XZ+EE(:,randperm(m-p)));

            [A,SIG] = tsdata_to_var(Uj,p);         % full regression
            LOGSIG = log(diag(SIG));
            [G,res] = var_to_autocov(A,SIG,qmax,acdecayfac);
            if res.error
                %fprintf(' *** bad VAR\n');
                nbad = nbad + 1;
                continue;
            end

            [~,SIGj] = autocov_to_var(G(jo,jo,:)); % reduced regression
            LOGSIGj = log(diag(SIGj));
            F(t,i,j) = LOGSIGj(ii)-LOGSIG(i);
            %fprintf('\n');

        end
        fprintf('done: %d of %d bad VARs\n',nbad,nperms);
    end

end
