function [G,itrs,relerr] = dlyap_itr(A,SIG,maxitrs,tol)

% [G,itrs,relerr] = dlyap_itr(A,SIG,maxitrs,tol)
%
% Solve discrete-time Lyapunov equation
%
%     G = A*G*A' + SIG
%
% (like 'dlyap', if you have the Control System Toolbox) using a simple
% iterative method. 'maxitrs' is the maximum iterations, 'tol' the
% tolerance (maximum relative error). The convergence rate depends on the
% dimension and maximum absolute eigenvalue of A (see e.g.
% utils/var_check.m).
%
% On return the caller should check for itrs == 0 (exceeded maximum
% iterations) and test relerr for non-convergence.

if nargin < 3 || isempty(maxitrs), maxitrs = 2000; end
if nargin < 4 || isempty(tol),     tol     = 1e-6; end

AT = A';
G = SIG;
for itrs = 1:maxitrs+1
    GG = G;
    G = A*GG*AT + SIG;
    relerr = norm(G-GG,'fro')/norm(G,'fro'); % Frobenius norm is much faster than the 2-norm
    if relerr < tol, break; end
end
if itrs > maxitrs, itrs = 0; end % exceeded maximum iterations
