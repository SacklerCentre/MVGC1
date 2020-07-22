function [G,res] = var_to_autocov_alt(A,SIG,acmaxlags,acdecayfac)

% [G,status,res.rho,res.aclags,res.acminlags] = var_to_autocov(A,SIG,acmaxlags,acdecayfac,maxlagfac)
%
% Return autocovariance sequence G for a VAR model with known coefficients A and
% residual covariance matrix SIG. For an n variable VAR(p), A must be n x n x p
% (i.e. p blocks of n x n) and SIG must be n x n symmetric positive-definite.
% This routine "reverse-solves" the Yule-Walker equations; it calls 'var_check'
% first to solve the associated 1-lag problem, essentially a discrete time
% Lyapunov equation.
%
% Various diagnostics and errors encountered in 'var_check' are passed through
% in the 'res' struct; in particular, res.rho conatins the spectral radius of
% the process. The VAR is stable iff res.rho < 1 (this may be considered a unit
% root test for stationarity). The res.error field MUST be checked by the
% caller: 0 signifies success, > 0 signifies an error with corresponding message
% in res.errmsg:
%
%  res.error   res.errmsg
%  --------------------------------------------------------------------
%     0        success
%     1        unstable VAR (has unit root)
%     2        supplied residuals covariance matrix not positive-definite
%     3        Lyapunov solver routine failed for some reason
%     4        calculated covariance matrix is not positive-definite
%
% The accuracy of the Lyapunov solution may be checked by res.actol, which
% should be very small.
%
% Further diagnostics are added to the res struct. For a stable VAR the
% autocovariance sequence decays approximately exponentially, by a factor equal
% to  res.rho. The minimum number of lags required to achieve the specified
% decay factor acdecayfac is calculated as 'res.acminlags'. The actual number of
% lags 'res.aclags' to which autocovariance is actually calculated is then set
% to the minimum of res.acminlags and the specified maximum number of lags,
% acmaxlags (if acmaxlags is not supplied - the recommended option - it defaults
% to res.acminlags).
%
% The caller should check that res.aclags >= res.acminlags; if this condition is
% not met, there is no guarantee that MVGCs - particularly in the spectral
% domain - will be accurate. However, if the spectral radius of the VAR model is
% close to one, so that res.acminlags is unfeasibly large, there may be no
% alternative [note that most Granger causality libraries effectively set
% res.aclags to the model order]. A typical check might be:
%
% if res.aclags < res.acminlags
%     fprintf(2,'WARNING: minimum %d lags required (decay factor = %e)\n',res.acminlags,realpow(res.rho,res.aclags));
% end

if nargin < 3 || isempty(acmaxlags), acmaxlags = 0; end % calculate maximum lags automatically (default)

if nargin < 4 || isempty(acdecayfac)
    acdecayfac = 1e-8; % autocovariance decay factor default
end

% check the VAR (don't need the 1-lag autocovariance - we're not going to
% use it).

res = var_check(A,SIG);

if res.error % something went wrong with the associated 1-lag calulation - bail out
    G = [];
    return;
end

% estimate required number of lags

p  = size(A,3);

res.acminlags = ceil(log(acdecayfac)/log(res.rho)); % minimum lags to achieve specified tolerance
if     acmaxlags < 0  % use exactly -acmaxlags lags (not encouraged, hence undocumented!)
    res.aclags = -acmaxlags;
elseif acmaxlags > 0  % use at most acmaxlags lags
    res.aclags = min(res.acminlags,acmaxlags);
else % acmaxlags == 0
    res.aclags = res.acminlags;
end
assert(res.aclags >= p,'number of lags is too small'); % lags must be at least number of VAR lags
% should test for res.aclags >= res.acminlags

q1 = res.aclags+1;

% calculate spectrum

G = cpsd_to_autocov(var_to_cpsd(A,SIG,q1));
