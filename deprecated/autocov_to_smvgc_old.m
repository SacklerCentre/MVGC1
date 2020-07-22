%% autocov_to_smvgc
%
% Calculate conditional frequency-domain MVGC (spectral multivariate Granger causality)
%
%% Syntax
%
%     [f,fres] = autocov_to_smvgc(G,x,y,fres)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     G          autocovariance sequence
%     x          vector of indices of target (causee) multi-variable
%     y          vector of indices of source (causal) multi-variable
%     fres       frequency resolution
%
% _output_
%
%     f          spectral Granger causality
%
%% Description
%
% Returns the frequency-domain (spectral) MVGC
%
% <<eq_smvgc.png>>
%
% from the variable |Y| (specified by the vector of indices |y|) to the
% variable |X| (specified by the vector of indices |x|), conditional on all
% other variables |Z| represented in |G|, for a stationary VAR process with
% autocovariance sequence G.
%
% Spectral causality is calculated up to the Nyqvist frequency at a
% resolution |fres|. If |fres| is not supplied it is calculated optimally
% as the number of autocovariance lags. Call |freqs =
% <sfreqs.html sfreqs>(fres,fs)|, where |fs| is the sampling
% rate, to get a corresponding vector |freqs| of frequencies on |[0,fs/2]|.
%
% In the conditional case, the algorithm works by transforming the
% autocovariance sequence for the full regression (see
% <autocov_xform.html |autocov_xform|>) to an autocovariance
% sequence for new |X,Z| variables defined as residuals of the reduced
% regression; thus a separate estimation step for the reduced regression,
% which is known to be problematic [2,*], is unnecessary, resulting in high
% efficiency and accuracy. See [1] for details.
%
%% References
%
% [1] L. Barnett and A. K. Seth, <matlab:openpdf('mvgc_doc.pdf') The MVGC
% Multivariate Granger Causality Toolbox>, _in preparation_, Aug. 2012.
%
% [2] Y. Chen, S. L. Bressler and M. Ding, "Frequency decomposition of
% conditional Granger causality and application to multivariate neural
% field potential data", _J. Neurosci. Methods_, 150, 2006. 
%
% [*] In our experience the "partition matrix" method in ref. [2] appears to be
% unsound, producing inaccurate results; hence we do not use it here.
%
%% See also
%
% <autocov_to_var.html |autocov_to_var|> |
% <var_to_cpsd.html |var_to_cpsd|> |
% <autocov_xform.html |autocov_xform|> |
% <var2trfun.html |var2trfun|> |
% <sfreqs.html |sfreqs|>
%
%% Copyright notice
%
% _(C) Lionel Barnett, 2012. See file <matlab:edit('license.txt') license.txt> in root directory
% for licensing terms._
%
%%

function [f,fres] = autocov_to_smvgc_old(G,x,y,fres)

[n,~,q1] = size(G);

x = x(:)'; % vectorise
y = y(:)'; % vectorise

assert(all(x >=1 & x <= n),'some x indices out of range');
assert(all(y >=1 & y <= n),'some y indices out of range');
assert(isempty(intersect(x,y)),'x and y indices must be distinct');

z = 1:n; z([x y]) = []; % indices of other variables (to condition out)

if nargin < 4 || isempty(fres);
    fres = q1;
end
h = fres+1;

f = zeros(1,h);

if isempty(z) % unconditional

    xy = [x y];
    G = G(xy,xy,:); % extract variables, rearrange
    nx = length(x);
    x = 1:nx;
    y = nx+1:n;

    [A,SIG] = autocov_to_var(G);                  % full regression
    [S,H] = var_to_cpsd(A,SIG,fres);              % spectrum & transfer function
    S = S(x,x,:);
    H = H(x,y,:);
    PSIG = SIG(y,y)-(SIG(y,x)/SIG(x,x))*SIG(x,y); % partial covariance
    for k = 1:h
        f(k) = log(real(det(S(:,:,k)))) - log(real(det(S(:,:,k)-H(:,:,k)*PSIG*H(:,:,k)')));
    end

else % conditional - transform autocov by reduced regression of x,z

    xz = [x z];
    [GR,~,SIGR] = autocov_xform(G,xz);
    nx = length(x);
    x = 1:nx;
    zy = nx+1:n;

    SRX = SIGR(x,x); % since reduced residuals serially uncorrelated, spectral density = covariance (flat spectrum) - and we only need x part

    % now do unconditional with transformed autocov

    [A,SIG] = autocov_to_var(GR);                     % transformed regression
    H = var2trfun(A,fres);                            % transfer function
    H = H(x,zy,:);
    PSIG = SIG(zy,zy)-(SIG(zy,x)/SIG(x,x))*SIG(x,zy); % partial covariance
    LDSRX = log(det(SRX));
    for k = 1:h
        f(k) = LDSRX-log(real(det(SRX-H(:,:,k)*PSIG*H(:,:,k)')));
    end

end
