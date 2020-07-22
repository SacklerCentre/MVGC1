%% var_to_tsdata
%
% Generate random multi-trial Gaussian VAR time series
%
%% Syntax
%
%     [X,E]        = var_to_tsdata(A,SIG,m,N,mtrunc)
%     [X,E,mtrunc] = var_to_tsdata(A,SIG,m,N,rho,decayfac)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     A          VAR coefficients matrix
%     SIG        residuals covariance matrix
%     m          number of observations per trial
%     N          number of trials (default: 1)
%     mtrunc     number of initial time steps to truncate
%     rho        VAR spectral radius (default: spectral radius of A)
%     decayfac   initial transients decay factor (default: 100)
%
% _output_
%
%     X          multi-trial Gaussian VAR time series
%     E          residuals time series
%     mtrunc     number of initial time steps truncated
%
%% Description
%
% Return |N| independent VAR time series of length |m| with coefficients |A| and
% iid Gaussian residuals |E|, with residuals covariance |SIG|:
%
% <<eq_var.png>>
%
% (where  [[ii_Sigma.png]] = |SIG|).
%
% _First form:_ if |mtrunc == 0| or |mtrunc >= 1| it is taken to be the the
% number of initial time steps to truncate.
%
% _Second form:_ if |0 < rho < 1| it is assumed to be the spectral radius of
% |A|; if omitted or empty it is set to the spectral radius of |A| (see function
% <var_specrad.html |var_specrad|>). The number of steps to assumed stationarity
% is then (over)estimated automatically and initial time steps truncated
% accordingly; set |decayfac| larger for longer settle time.
%
%
%% References
%
% [1] L. Barnett and A. K. Seth, <matlab:openpdf('mvgc_doc.pdf') The MVGC
% Multivariate Granger Causality Toolbox>, _in preparation_, Aug. 2012.
%
%% See also
%
% <var_specrad.html |var_specrad|> |
% <genvar.html |genvar|>
%
%% Copyright notice
%
% _(C) Lionel Barnett, 2012. See file <matlab:edit('license.txt')
% license.txt> in root directory for licensing terms._
%
%%

function [X,E,mtrunc] = var_to_tsdata_old(A,SIG,m,N,mtrho,decayfac)

if nargin < 4 || isempty(N), N = 1; end % single trial

if nargin < 5 || isempty(mtrho)
    mtrho = var_specrad(A);
    assert(mtrho < 1,'unstable VAR');
end

if nargin < 6 || isempty(decayfac)
    decayfac = 100; % should be more than enough...
end

[C,cholp] = chol(SIG,'lower');
assert(cholp == 0,'covariance matrix not positive-definite');

n = size(A,1);

if mtrho > 0 && mtrho < 1 % assume mtrho is the spectral radius of A
    mtrunc = round((log(eps)-decayfac)/log(mtrho)); % enough time for autocovariance to decay to fp accuracy (and then some)
else                      % assume mtrho is the number of  time steps to truncate
    mtrunc = round(mtrho);
end

if N > 1 % multi-trial

    X = zeros(n,m,N);
    if nargout > 1
        E = zeros(n,m,N);
        for r = 1:N
            [X(:,:,r),E(:,:,r)] = genvar(A,C*randn(n,m+mtrunc),mtrunc);
        end
    else
        for r = 1:N
            X(:,:,r) = genvar(A,C*randn(n,m+mtrunc),mtrunc);
        end
    end

else

    if nargout > 1
        [X,E] = genvar(A,C*randn(n,m+mtrunc),mtrunc);
    else
        X = genvar(A,C*randn(n,m+mtrunc),mtrunc);
    end

end
