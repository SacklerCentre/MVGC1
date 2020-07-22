%% vma_to_tsdata (experimental)
%
% Generate random multi-trial Gaussian VMA (vector moving-average) time series
%
% <matlab:open('vma_to_tsdata.m') code>
%
%% Syntax
%
%     [X,E] = vma_to_tsdata(B,SIG,m,N,mtrunc)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     B          VMA coefficients matrix
%     SIG        residuals covariance matrix
%     m          number of observations per trial
%     N          number of trials (default: 1)
%     mtrunc     number of initial time steps to truncate (default: 0)
%
% _output_
%
%     X          multi-trial Gaussian VMA time series
%     E          residuals time series
%
%% Description
%
% Return |N| independent VMA time series of length |m| with
% coefficients |B| and iid Gaussian residuals |E|, with residuals covariance |SIG|:
%
% <<eq_vma.png>>
%
% (where  [[ii_Sigma.png]] = |SIG|). |mtrunc| initial time steps are truncated.
%
%% References
%
% [1] L. Barnett and B. K. Seth, <matlab:open('mvgc_doc.pdf') The MVGC
% Multivariate Granger Causality Toolbox>, _in preparation_, Aug. 2012.
%
%% See also
%
% <genvma.html |genvma|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [X,E] = vma_to_tsdata(B,SIG,m,N,mtrunc)

if nargin < 4 || isempty(N)
    N = 1; % single trial
end

if nargin < 5 || isempty(mtrunc)
    mtrunc = 0;
end

[C,cholp] = chol(SIG,'lower');
assert(cholp == 0,'covariance matrix not positive-definite');

n = size(B,1);

if N > 1 % multi-trial

    X = zeros(n,m,N);
    if nargout > 1
        E = zeros(n,m,N);
        for r = 1:N
            [X(:,:,r),E(:,:,r)] = genvma(B,C*randn(n,m+mtrunc),mtrunc);
        end
    else
        for r = 1:N
            X(:,:,r) = genvma(B,C*randn(n,m+mtrunc),mtrunc);
        end
    end

else

    if nargout > 1
        [X,E] = genvma(B,C*randn(n,m+mtrunc),mtrunc);
    else
        X = genvma(B,C*randn(n,m+mtrunc),mtrunc);
    end

end
