%% tsdata_to_var
%
% Fit VAR model to multi-trial, multivariate time series data
%
% <matlab:open('tsdata_to_var.m') code>
%
%% Syntax
%
%     [A,SIG,E] = tsdata_to_var(X,p,regmode)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     X          multi-trial time series data
%     p          model order (number of lags)
%     regmode    regression mode: 'LWR' (default) or 'OLS'
%
% _output_
%
%     A          VAR coefficients matrix
%     SIG        residuals covariance matrix
%     E          residuals time series
%
%% Description
%
% Returns VAR coefficients |A| and (optionally) residuals covariance matrix
% |SIG| and serially uncorrelated residuals |E| for the |p|-lag autoregression
%
% <<eq_var.png>>
%
% (where  [[ii_Sigma.png]] = |SIG|) of a stationary multivariate process
% |X|. |X| may contain single- or multi-trial multivariate time series
% data. The regression mode is set by the |regmode| parameter, which may be
% |'LWR'| (default) or |'OLS'|. The former uses Morf's version of the LWR
% algorithm [1,2] while the latter calculates the OLS solution to the
% regression via QR decomposition.
%
% *_Note_*: If the regressions are rank-deficient or ill-conditioned then A may
% be "bad" (i.e. will contain a |NaN| or |Inf|; see <isbad.html |isbad|>) and/or
% warnings will may be issued. The caller should test for both these
% possibilities, by calls to <isbad.html |isbad|> and <warn_supp.html
% |warn_supp|> ... <warn_test.html |warn_test|> respectively. Possible causes
% are non-stationarity and/or colinearity in the data.
%
% The caller should also, at the very least, check the _spectral radius_ of the
% returned VAR coefficients (see <var_specrad |var_specrad|>) to ensure that the
% coefficients define a stable VAR [1]. (This is calculated, along with other
% relevant information, in the routine <var_to_autocov.html |var_to_autocov|>,
% which will typically be called subsequent to this function, and may be tested
% by a call to <var_acinfo.html |var_acinfo|>).
%
% *_Note_*:  LWR implementation patched March 2018 to correct an initialisation
% error (many thanks to Gonzalo Camba-Mendez for the heads-up). The error meant
% that the algorithm was not strictly "minimum phase", so that coefficien
% estimates were potentially unstable. (Note that the OLS algorithm, which
% returns (pseudo) maximum-likelihood estimates, is not guaranteed to return
% stable estimates.)
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% [2] M. Morf, A. Viera, D. T. L. Lee and T. Kailath, "Recursive Multichannel
% Maximum Entropy Spectral Estimation", _IEEE Trans. Geosci. Elec._, 16(2), 1978.
%
%% See also
%
% <var_specrad.html |var_specrad|> |
% <var_to_autocov.html |var_to_autocov|> |
% <warn_supp.html |warn_supp|> |
% <warn_test.html |warn_test|> |
% <isbad.html |isbad|> |
% <var_acinfo.html |var_acinfo|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [A,SIG,E] = tsdata_to_var(X,p,regmode)

if nargin < 3 || isempty(regmode), regmode = 'LWR'; end

[n,m,N] = size(X);
assert(p < m,'too many lags');

p1 = p+1;
pn = p*n;
p1n = p1*n;

A   = NaN; % assure a "bad" return value if anything goes wrong (see routine 'isbad')
SIG = NaN; % assure a "bad" return value if anything goes wrong (see routine 'isbad')
E   = NaN; % assure a "bad" return value if anything goes wrong (see routine 'isbad')

X = demean(X); % no constant term

if  strcmpi(regmode,'OLS') % OLS (QR decomposition)

    M = N*(m-p);

    % stack lags

    X0 = reshape(X(:,p1:m,:),n,M); % concatenate trials for unlagged observations
    XL = zeros(n,p,M);
    for k = 1:p
        XL(:,k,:) = reshape(X(:,p1-k:m-k,:),n,M); % concatenate trials for k-lagged observations
    end
    XL = reshape(XL,pn,M);         % stack lags

    A = X0/XL;                     % OLS using QR decomposition
    if isbad(A); return; end       % something went badly wrong

    if nargout > 1
        E   = X0-A*XL;             % residuals
        SIG = (E*E')/(M-1);        % residuals covariance matrix
        E   = reshape(E,n,m-p,N);  % put residuals back into per-trial form
    end

    A = reshape(A,n,n,p);          % so A(:,:,k) is the k-lag coefficients matrix

elseif strcmpi(regmode,'LWR') % LWR (Morf): patched, March 2018 to correct initialisation error (many thanks to Gonzalo Camba-Mendez for the heads-up)

    I = eye(n);

    % store lags

    XX = zeros(n,p1,m+p,N);
    for k = 0:p
        XX(:,k+1,k+1:k+m,:) = X; % k-lagged observations
    end

    % initialise recursion

    EE = reshape(X,n,N*m);
    IC = inv(chol(EE*EE','lower')); % inverse covariance square root

    k  = 1;
    kn = k*n;
    M  = N*(m-k);
    kk = 1:k;
    kf = 1:kn;         % forward  indices
    kb = p1n-kn+1:p1n; % backward indices

    AF = zeros(n,p1n); AF(:,kf) = IC; % forward  AR coefficients
    AB = zeros(n,p1n); AB(:,kb) = IC; % backward AR coefficients (reversed compared with [2])

    % LWR recursion

    while k <= p

        EF = AF(:,kf)*reshape(XX(:,kk,k+1:m,:),kn,M); % forward  prediction errors
        EB = AB(:,kb)*reshape(XX(:,kk,k:m-1,:),kn,M); % backward prediction errors

        R = (chol(EF*EF','lower')\EF)*(chol(EB*EB','lower')\EB)'; % normalised reflection coefficients

        k  = k+1;
        kn = k*n;
        M  = N*(m-k);
		kk = 1:k;
        kf = 1:kn;
        kb = p1n-kn+1:p1n;

        AFPREV = AF(:,kf);
        ABPREV = AB(:,kb);

        AF(:,kf) = chol(I-R*R','lower')\(AFPREV-R*ABPREV);
        AB(:,kb) = chol(I-R'*R,'lower')\(ABPREV-R'*AFPREV);

    end

	A0 = AF(:,1:n);
    A = reshape(-A0\AF(:,n+1:p1n),n,n,p);
    if isbad(A); return; end       % something went badly wrong

    if nargout > 1
		M   = N*(m-p);      % residuals lose p lags
		E   = A0\EF;        % residuals
        SIG = (E*E')/(M-1); % residuals covariance matrix (unbiased estimator)
        if nargout > 2      % align residuals per-trial with data (lose p lags)
			E = cat(2,nan(n,p,N),reshape(E,n,m-p,N));
		end
    end

else
    error('bad regression mode ''%s''',regmode);
end
