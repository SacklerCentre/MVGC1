function err = var_check(A,SIG)

% err = var_check(A,SIG)
%
% Calculates the covariance matrix for the VAR(1) associated with the p-lag VAR
% specified by coefficients matrix A and residuals covariance matrix SIG. For an
% n variable VAR(p), A must be n x n x p  (i.e. p blocks of n x n) and SIG must
% be n x n symmetric positive definite. The returned covariance matrix G1 is p*n
% x p*n, organised in blocks of n x n, so that the i,j-th block corresponds to
% j-i lags (the top n rows are also the autocovariance sequence up to p-1 lags
% for the p-lag VAR specified by A, SIG; see 'var_to_autocov').
%
% Various results and diagnostics are returned in the 'res' struct. In
% particular, the associated 1-lag problem is checked; res.rho returns its
% spectral radius - the VAR is stable iff res.rho < 1 (this may be
% considered a unit root test for stationarity). The 1-lag autocovariance
% matrix is also checked for positive-definitiveness. If this check fails,
% it may be an indication of an ill-conditioned VAR (possibly because
% residuals variances are too small, and/or the process is borderline
% stationary). In short, the res.error field must be checked by the caller:
% 0 signifies success, > 0 signifies an error with corresponding message in
% res.errmsg:
%
%  res.error   res.errmsg
%  --------------------------------------------------------------------
%     0        -
%     1        unstable VAR (has unit root)
%     2        supplied residuals covariance matrix not positive-definite
%     3        Lyapunov solver routine failed for some reason
%     4        1-lag covariance matrix is not positive-definite
%
% The accuracy of the Lyapunov solution may be checked by res.actol, which
% should be very small.
%
% If the Matlab function 'dlyap' (Control System Toolbox) is not available, a
% slower routine based on Schur decomposition is used; see utils/dlyap_schur.m.
% [NOTE: In future releases we hope to replace this function with a SLICOT
% routine, licensing issues permitting.] In fact, the slower routine may
% sometimes be more stable, and may be selected by setting the 'schur_dlyap'
% flag to true.

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
pn1 = (p-1)*n;

[nn1,nn2] = size(SIG);
assert(nn1 == nn2,'residuals covariance matrix not square');
assert(nn1 == n  ,'residuals covariance matrix doesn''t match VAR coefficients matrix');

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

res = rmfield(res,'errmsg'); % remove error message field if no errors
