%% parcov
%
% Calculate partial covariance
%
% <matlab:open('parcov.m') code>
%
%% Syntax
%
%     P = parcov(V,x,y)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     V        covariance matrix
%     x        vector of indices of target multi-variable
%     y        vector of indices of conditioning multi-variable
%
% _output_
%
%     P        partial covariance matrix
%
%% Description
%
% Given a (symmetric, positive-definite) covariance matrix [[ii_Sigma.png]],
% calculates the (symmetric, positive-definite) partial covariance
%
% <<eq_parcov.png>>
%
% *_WARNING_*: this function does no parameter or error checking, as it will
% generally be called in computationally intensive loops, in particular for
% spectral GC calculation (<autocov_to_smvgc.html |autocov_to_smvgc|>,
% <autocov_to_spwcgc.html |autocov_to_spwcgc|>). _The caller must ensure that
% parameters are in order and that_ [[ii_Sigma.png]] _is positive-definite_!
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
%% See also
%
% <autocov_to_smvgc.html |autocov_to_smvgc|> |
% <autocov_to_spwcgc.html |autocov_to_spwcgc|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function P = parcov(V,x,y)

if isscalar(y)
    U = sqrt(V(y,y))\V(y,x);   % faster if y 1-dim
else
    U = linsolve(chol(V(y,y),'lower'),V(y,x),struct('LT',true)); % 'linsolve' is more efficient than 'rdivide' (we know the Cholesky factor is lower-triangular)
end
P = V(x,x)-U'*U;
