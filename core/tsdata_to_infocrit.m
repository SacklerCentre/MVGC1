%% tsdata_to_infocrit
%
% Calculate Akaike and Bayesian information criteria for VAR models from time series
% data
%
% <matlab:open('tsdata_to_infocrit.m') code>
%
%% Syntax
%
%     [aic,bic,moaic,mobic] = tsdata_to_infocrit(X,morder,regmode,verb)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
% time sries data
% _input_
%
%     X          multi-trial time series data
%     morder     maximum VAR model order or vector of model orders (see Description)
%     regmode    regression mode: 'LWR' (default) or 'OLS'
%     verb       verbosity flag: true (default) or false
%
% _output_
%
%     aic        vector of AIC values
%     bic        vector of BIC values
%     moaic      optimal model order according to AIC
%     mobic      optimal model order according to BIC
%
%% Description
%
% Calculates Akaike (AIC) and Bayesian (BIC) information criteria for VAR models
% from time series data |X|, which may be single- or multi-trial. Information
% criteria are returned in the (column) vectors |aic|, |bic|; see <infocrit.html
% |infocrit|>. Optimal model orders according to the respective criteria are
% returned in |moaic|, |mobic|.
%
% The regression mode is set by the |regmode| parameter, which may be |'LWR'|
% (default) or |'OLS'|. The former uses Morf's version of the LWR algorithm
% [1,2] while the latter calculates the OLS solution to regressions via QR
% decomposition. The LWR algorithm should be faster, as regressions are
% performed recursively so that the information criteria may be calculated "on
% the fly".
%
% The LWR algorithm must fit a VAR model at all model orders up to the maximum;
% in this case the |morder| parameter must be a scalar representing the maximum
% model order. For the OLS mode, |morder| may again be a scalar representing the
% maximum model order, or it may be a vector of model orders.
%
% If the |verb| (verbosity) flag is set (default), progress is reported;
% otherwise only warnings are reported.
%
% If a VAR parameter estimate fails - because of ill-conditioned regression, or
% the residuals covariance matrix comes out as non positive-definite - a |NaN|
% is returned for that model order. This may be an indication that the time
% series data on which the autocovariance sequence estimate was based is not
% sufficiently stationary, is not of sufficient length, has (possibly lagged)
% colinearities, or has a highly skewed distribution.
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
% Maximum Entropy Spectral Estimation", _IEEE Trans. Geosci. Elec._, 16(2),
% 1978.
%
%% See also
%
% <infocrit.html |infocrit|> |
% <tsdata_to_var.html |tsdata_to_var|> |
% <mvgc_demo.html |mvgc_demo|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [aic,bic,moaic,mobic] = tsdata_to_infocrit(X,morder,regmode,verb)

if nargin < 3 || isempty(regmode), regmode = 'LWR'; end
if nargin < 4 || isempty(verb),    verb    = true;  end

[n,m,N] = size(X);

X = demean(X,2); % no constant term, normalise

% store lags

q = max(morder);
q1 = q+1;
XX = zeros(n,q1,m+q,N);
for k = 0:q
    XX(:,k+1,k+1:k+m,:) = X; % k-lagged observations
end

if  strcmpi(regmode,'OLS') % OLS (QR decomposition)

    assert(isvector(morder) || isscalar(morder) && all(morder > 0),'model orders must be a positive vector or scalar');
    if isscalar(morder), morder = 1:morder; end

    nummo = length(morder);

    aic = nan(nummo,1);
    bic = nan(nummo,1);

    % loop through model orders

    for i = 1:nummo

        q = morder(i);

        if verb, fprintf('model order = %d',q); end

        if q >= m
            if ~verb, fprintf('model order = %d',q); end
            fprintf(2,'  WARNING: model order too large (must be < %d)\n',m);
            continue
        end

        q1 = q+1;
        M = N*(m-q);
        nq = n*q;

        % stack lags

        X0 = reshape(XX(:,1,q1:m,:),n,M);
        XL = reshape(XX(:,2:q1,q1:m,:),nq,M);

        owstate = warn_supp;
        A = X0/XL;                                         % OLS using QR decomposition
        [iw,~,wid]= warn_test(owstate,[],false);
        if isbad(A)                                        % something went badly wrong
            if ~verb, fprintf('model order = %d',q); end
            fprintf(2,'  WARNING: VAR estimation failed\n');
            continue % show-stopper
        end
        if iw % rank-deficient?
            if ~verb, fprintf('model order = %d',q); end
            fprintf(2,'  WARNING: VAR estimation may be problematic (%s)',wid);
            % not necessarily a show-stopper - carry on
        end

        E   = X0-A*XL;            % residuals
        DSIG = det((E*E')/(M-1)); % residuals covariance matrix determinant
        if DSIG <= 0
            if ~iw, if ~verb, fprintf('model order = %d',q); end; end
            fprintf(2,'  WARNING: residuals covariance not positive definite\n');
            continue % show-stopper
        end

        [aic(i),bic(i)] = infocrit(-(M/2)*log(DSIG),q*n*n,M); % -(M/2)*log(DSIG) is max log-likelihood
        if verb, fprintf(1,'\n'); end
    end

elseif strcmpi(regmode,'LWR') % LWR (Morf): patched, March 2018 to correct initialisation error (many thanks to Gonzalo Camba-Mendez for the heads-up)

    assert(isscalar(morder) && morder > 0,'model order must be a positive scalar');

    aic = nan(morder,1);
    bic = nan(morder,1);

	q = morder;
	q1 = q+1;
	qn = q*n;
	q1n = q1*n;

	I = eye(n);

    % initialise recursion

	M = N*m;
	E = reshape(X,n,M);
    IC = inv(chol(E*E','lower')); % inverse covariance square root

    k  = 1;
    kn = k*n;
    M  = N*(m-k);
    kk = 1:k;
    kf = 1:kn;         % forward  indices
    kb = q1n-kn+1:q1n; % backward indices

    AF = zeros(n,q1n); AF(:,kf) = IC; % forward  AR coefficients
    AB = zeros(n,q1n); AB(:,kb) = IC; % backward AR coefficients (reversed compared with [2])

    % LWR recursion

    while k <= q

        if verb, fprintf('model order = %d',k); end

        EF = AF(:,kf)*reshape(XX(:,kk,k+1:m,:),kn,M); % forward  prediction errors
        EB = AB(:,kb)*reshape(XX(:,kk,k:m-1,:),kn,M); % backward prediction errors

        wstate = warning('off','all'); lastwarn('');
        R = (chol(EF*EF','lower')\EF)*(chol(EB*EB','lower')\EB)'; % normalised reflection coefficients
        wmsg = lastwarn; warning(wstate);
        if ~isempty(wmsg) % rank-deficient?
            if ~verb, fprintf('model order = %d',k); end
            fprintf(2,'  WARNING: reflection coefficients estimation may be problematic (%s)',wmsg);
            % not necessarily a show-stopper - carry on
        end
        if isbad(R) % something went badly wrong
            if ~verb, fprintf('model order = %d',k); end
            fprintf(2,'  WARNING: reflection coefficients estimation failed\n');
            continue % show-stopper
        end

		kp = k;
        k  = k+1;
        kn = k*n;
        M  = N*(m-k);
		kk = 1:k;
        kf = 1:kn;
        kb = q1n-kn+1:q1n;

        AFPREV = AF(:,kf);
        ABPREV = AB(:,kb);

        wstate = warning('off','all'); lastwarn('');
        AF(:,kf) = chol(I-R*R','lower')\(AFPREV-R*ABPREV);
        AB(:,kb) = chol(I-R'*R,'lower')\(ABPREV-R'*AFPREV);
        wmsg = lastwarn; warning(wstate);
        if ~isempty(wmsg) % rank-deficient?
            if ~verb, fprintf('model order = %d',kp); end
            fprintf(2,'  WARNING: forward/backward VAR coefficients estimation may be problematic (%s)',wmsg);
            % not necessarily a show-stopper - carry on
        end
        if isbad(AF) || isbad(AB) % something went badly wrong
            if ~verb, fprintf('model order = %d',kp); end
            fprintf(2,'  WARNING: forward/backward VAR coefficients estimation failed\n');
            continue % show-stopper
        end

        E = AF(:,1:n)\EF; % residuals

        DSIG = det((E*E')/(M-1)); % residuals covariance matrix determinant
        if DSIG <= 0
            if ~verb, fprintf('model order = %d',kp); end
            fprintf(2,'  WARNING: residuals covariance not positive definite\n');
            continue % show-stopper
        end

        [aic(kp),bic(kp)] = infocrit(-(M/2)*log(DSIG),kp*n*n,M); % -(M/2)*log(DSIG) is max log-likelihood
        if verb, fprintf(1,'\n'); end
    end

else
    error('bad regression mode ''%s''',regmode);
end

if nargout > 2 % calculate optimal model orders (NaNs are ignored)
    if isscalar(morder), morder = 1:morder; end
    [~,idx] = min(aic); moaic = morder(idx);
    [~,idx] = min(bic); mobic = morder(idx);
end
