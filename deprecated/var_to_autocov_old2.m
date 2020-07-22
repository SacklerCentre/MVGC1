%% var_to_autocov
%
% Return autocovariance sequence for a VAR model
%
% <matlab:open('var_to_autocov.m') code>
%
%% Syntax
%
%     [G,res] = var_to_autocov(A,SIG,acmaxlags,acdectol,dlyap_alg,maxiters,maxrelerr)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     A          VAR coefficients matrix
%     SIG        residuals covariance matrix
%     acmaxlags  maximum autocovariance lags to calculate or (default) zero for automatic calculation
%     acdectol   autocovariance decay tolerance (default: 1e-8)
%     dlyap_alg  Lyapunov solver algorithm: 'builtin' (default), 'Schur' or 'aitr'; see Description for meanings)
%     maxiters   maximum iterations if using iterative algorithm; see 'dlyap_aitr' for defaults
%     maxrelerr  maximum relative error if using iterative algorithm; see 'dlyap_aitr' for defaults
%
% _output_
%
%     G          autocovariance sequence
%     res        results structure, with fields (some may not be present):
%         error      error number (0 for no error)
%         errmsg     error message string
%         rho        VAR spectral radius
%         iterations number of iterations performed if using iterative algorithm
%         acrelerr   relative error of associated 1-lag solution
%         acminlags  minimum lags required to achieve specified autocovariance decay factor
%         aclags     actual number of autocovariance lags calculated
%         WARNINGn   warnings
%
%% Description
%
% Returns autocovariance sequence |G| defined as [[ii_acseq.png]]
% for a VAR model with coefficients |A| and (positive-definite) residual
% covariance matrix |SIG|, by "reverse-solving" the Yule-Walker equations
%
% <<eq_yweqs.png>>
%
% (where  [[ii_Sigma.png]] = |SIG|). The algorithm solves the associated
% 1-lag problem - a discrete-time Lyapunov equation - and then calculates
% higher lags recursively [1].
%
% Errors and diagnostics are returned in the |res| struct. Possible errors are
%
%     res.error      res.errmsg
%     ----------------------------------------------------------------
%         0          (no error, no message)
%         1          unstable VAR (has unit root)
%         2          residuals covariance matrix not positive-definite
%         3          Lyapunov equation solver failed for some reason
%         4          1-lag covariance matrix not positive-definite
%     ----------------------------------------------------------------
%
% For a stable VAR the the spectral radius |res.rho|  (see
% <var_specrad.html |var_specrad|>) must be < 1; this may be
% considered a unit root test for stationarity [1]. Then the autocovariance
% sequence decays approximately exponentially, by a factor equal to
% |res.rho|. The minimum number of lags required to achieve the specified
% autocovariance decay tolerance |acdectol| is calculated as
% |res.acminlags|, so that |res.rho^res.acminlags < acdectol|. The actual
% number of lags |res.aclags| to which autocovariance is calculated is then
% set to the minimum of |res.acminlags| and the specified maximum number of
% lags, |acmaxlags| (if |acmaxlags| is not supplied - the recommended
% option - it defaults to |res.acminlags|). A warning is issued if
% |res.aclags < res.acminlags|. In this case there is no guarantee that
% MVGCs - particularly in the spectral domain - will be accurate. However,
% if the spectral radius of the VAR model is close to 1, so that
% |res.acminlags| is unfeasibly large, there may be no alternative [note
% that most Granger causality libraries effectively set |res.aclags| to the
% model order]. The calculated 1-lag autocovariance matrix is also checked
% for positive-definitiveness. If this check fails, it may be an indication
% of an ill-conditioned VAR (possibly because residuals variances are too
% small, and/or the process is borderline stationary). In short, the
% |res.error| field _must_ be checked by the caller: 0 signifies success, >
% 0 signifies an error with corresponding message in |res.errmsg|.
% Attention should also be paid to any warnings in |res.WARNINGn|.
%
% The string |dlyap_alg| specifies the Lyapunov solver algorithm. It may be set
% to |'builtin'| for the Control System Toolbox <matlab:doc('dlyap') |dlyap|>
% solver routine (if available), |'Schur'| for a roughly equivalent (but slower)
% scripted algorithm based on Schur decomposition (see <dlyap_schur.html
% |dlyap_schur|>), or |'aitr'| for an experimental "accelerated" iterative
% algorithm (see <dlyap_aitr.html |dlyap_aitr|>). The current default is
% |'builtin'| if the Control System Toolbox |dlyap| is available, |'Schur'|
% otherwise. [ _Note:_ in future releases we hope to replace the |'Schur'|
% algorithm with the appropriate <http://www.slicot.org/ SLICOT> routine,
% licensing issues permitting. ]
%
%% References
%
% [1] L. Barnett and A. K. Seth, <matlab:open('mvgc_doc.pdf') The MVGC
% Multivariate Granger Causality Toolbox>, _in preparation_, Aug. 2012.
%
%% See also
%
% <var_specrad.html |var_specrad|> |
% <matlab:doc('dlyap') |dlyap|> |
% <dlyap_schur.html |dlyap_schur|> |
% <dlyap_aitr.html |dlyap_aitr|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [G,res] = var_to_autocov(A,SIG,acmaxlags,acdectol,dlyap_alg,maxiters,maxrelerr)

global have_dlyap;

% default parameters

if nargin < 3 || isempty(acmaxlags), acmaxlags = 0;    end % calculate maximum lags automatically
if nargin < 4 || isempty(acdectol),  acdectol  = 1e-8; end % autocovariance decay tolerance
if nargin < 5 || isempty(dlyap_alg)
    if have_dlyap % got Control System Toolbox dlyap
        dlyap_alg = 'builtin';
    else
        dlyap_alg = 'schur'; % see utils/dlyap_schur.m.
    end
end

% iterative algorithm only: ensure defaults for utils/dlyap_aitr.m.

if nargin < 6, maxiters  = []; end
if nargin < 7, maxrelerr = []; end

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
pn1 = (p-1)*n;

[nn1,nn2] = size(SIG);
assert(nn1 == nn2,'residuals covariance matrix not square');
assert(nn1 == n  ,'residuals covariance matrix doesn''t match VAR coefficients matrix');

res.error  = 0;
res.errmsg = '';
G          = [];

% construct VAR coefficients for 1-lag problem

A1 = [reshape(A,n,p*n); eye(pn1) zeros(pn1,n)];

% calculate spectral radius

res.rho = max(abs(eig(A1)));

if res.rho >= 1
    res.error = 1;
    res.errmsg = 'unstable VAR (unit root)';
    return
end

% construct residual covariances for 1-lag problem

if ~isposdef(SIG);
    res.error = 2;
    res.errmsg = 'residuals covariance matrix not positive-definite';
    return;
end

SIG1 = [SIG zeros(n,pn1); zeros(pn1,n) zeros(pn1)];

% solve the Lyapunov equation for the 1-lag covariance matrix

try
    switch lower(dlyap_alg)
        case 'builtin'
            assert(have_dlyap,'builtin Lyapunov solver (Control System Toolbox) not available');
            %G1 = dlyap(A1,SIG1);          % dlyap seems to work better here without balancing, which seems to break positive-definitiveness
            G1 = lyapslv('D',A1,[],-SIG1); % sometimes. However lyapslv is not an official interface, so this could conceivably break in future.
        case 'schur'
            G1 = dlyap_schur(A1,SIG1);     % slow (but pretty much same algorithm as builtin 'dlyap')
        case 'aitr'
            [G1,res.iterations] = dlyap_aitr(A1,SIG1,maxiters,maxrelerr); % experimental: fast, but needs more testing
        otherwise
            error('unknown Lyapunov solver ''%s''',dlyap_alg);
    end
catch except
    res.error = 3;
    res.errmsg = ['Lyapunov equation solver failed: ' except.message];
    return
end

res.acrelerr = norm(A1*G1*A1'-G1+SIG1)/norm(SIG1); % this should be small (see below)

% estimate number of autocov lags

res.acminlags = ceil(log(acdectol)/log(res.rho)); % minimum lags to achieve specified tolerance

if     acmaxlags < 0  % use exactly -acmaxlags lags (not encouraged, hence undocumented!)
    res.aclags = -acmaxlags;
elseif acmaxlags > 0  % use at most acmaxlags lags
    res.aclags = min(res.acminlags,acmaxlags);
else                  % acmaxlags == 0 - use minimum acceptable lags (recommended)
    res.aclags = res.acminlags;
end

maxacrelerr = 1e-8;
if res.acrelerr > maxacrelerr
    res.WARNING1 = 'large relative error for Lyapunov equation solution';
    if res.aclags < res.acminlags
        res.WARNING2 = 'too few autocovariance lags to achieve specified autocovariance decay tolerance';
    end
else
    if res.aclags < res.acminlags
        res.WARNING1 = 'too few autocovariance lags to achieve specified autocovariance decay tolerance';
    end
end

if ~isposdef(G1);
    res.error = 4;
    res.errmsg = '1-lag covariance matrix not positive-definite';
    return
end

res = rmfield(res,'errmsg'); % no errors if we got here, so remove error message field

q = res.aclags;
q1 = q+1;

% calculate recursively from 1-lag solution (which supplies up to p-1 lags), from p lags up to q

[n,~,p]  = size(A);
assert(res.aclags >= p,'number of lags is too small'); % lags must be at least number of VAR lags
pn = p*n;
G = cat(3,reshape(G1(1:n,:),n,n,p),zeros(n,n,q1-p));   % autocov forward  sequence
B = [zeros((q1-p)*n,n); G1(:,end-n+1:end)];            % autocov backward sequence
A = reshape(A,n,pn);                                   % coefficients
for k = p:q
    r = q1-k;
    G(:,:,k+1) = A*B(r*n+1:r*n+pn,:);
    B((r-1)*n+1:r*n,:) = G(:,:,k+1);
end
