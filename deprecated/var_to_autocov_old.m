function [G,res] = var_to_autocov(A,SIG,acmaxlags,acdecayfac,schur_dlyap)

% [G,res] = var_to_autocov(A,SIG,acmaxlags,acdecayfac,schur_dlyap)
%
% Return autocovariance sequence G for a VAR model with known coefficients
% A and residual covariance matrix SIG. For an n variable VAR(p), A must be
% n x n x p (i.e. p blocks of n x n) and SIG must be n x n symmetric
% positive-definite. This routine "reverse-solves" the Yule-Walker
% equations; it first calls 'var_check', which solves the associated 1-lag
% problem - a discrete time Lyapunov equation - then calculates higher lags
% recursively. (If the Lyapunov solver is not available a slower iterative
% method is used; see utils/var_check.m.)
%
% Various diagnostics and errors encountered in 'var_check' are passed through
% in the 'res' struct. The res.error field MUST be checked by the caller -
% ignoring an error will most likely cause failure of subsequent routines that
% attempt to use the returned autocovariance sequence; see 'var_check' for
% details.
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
%
% For the meaning of the 'schur_dlyap' flag, see utils/var_check.m.

% default parameters

if nargin < 3 || isempty(acmaxlags),   acmaxlags  = 0;    end % calculate maximum lags automatically
if nargin < 4 || isempty(acdecayfac),  acdecayfac = 1e-8; end % autocovariance decay factor
if nargin < 5,                         schur_dlyap = [];  end % to pick up var_check() default

% check the VAR/solve 1-lag problem

[res,G1] = var_check(A,SIG,schur_dlyap); % check VAR parameters and solve the 1-lag problem

if res.error % bad VAR - bail out
    G = [];
    return;
end

% estimate number of autocov lags

res.acminlags = ceil(log(acdecayfac)/log(res.rho)); % minimum lags to achieve specified tolerance

if     acmaxlags < 0  % use exactly -acmaxlags lags (not encouraged, hence undocumented!)
    res.aclags = -acmaxlags;
elseif acmaxlags > 0  % use at most acmaxlags lags
    res.aclags = min(res.acminlags,acmaxlags);
else                  % acmaxlags == 0 - use minimum acceptable lags (recommended)
    res.aclags = res.acminlags;
end

% NOTE: should test for res.aclags >= res.acminlags

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
