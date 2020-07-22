%% autocov_to_smvgc_pwc
%
% Calculate pairwise-conditional frequency-domain MVGCs (spectral multivariate Granger causalites)
%
%% Syntax
%
%     [f,fres] = autocov_to_smvgc_pwc(G,fres)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     G          autocovariance sequence
%     fres       frequency resolution
%
% _output_
%
%     F          spectral Granger causality matrix
%
%% Description
%
% Returns the  matrix |f| of pairwise-conditional frequency-domain (spectral) MVGCs
%
% <<eq_smvgc_pwc.png>>
%
% (where |[ij]| denotes omission of the |ij|-th variables) between all
% pairs of variables |i $\ne$ j| represented in |G|, for a stationary VAR
% process with autocovariance sequence |G|. The first index |i| of
% |f| is the target (causee) variable, the second |j| the source (causal)
% variable and the third indexes the frequency. See ref. [1] for details.
%
% Spectral causality is calculated up to the Nyqvist frequency at a
% resolution |fres|. If |fres| is not supplied it is calculated optimally
% as the number of autocovariance lags. Call |freqs =
% <sfreqs.html sfreqs>(fres,fs)|, where |fs| is the sampling
% rate, to get a corresponding vector |freqs| of frequencies on |[0,fs/2]|.
%
% For details of the algorithm, see <autocov_to_smvgc.html |autocov_to_smvgc|> and [1].
%
%% References
%
% [1] L. Barnett and A. K. Seth, <matlab:openpdf('mvgc_doc.pdf') The MVGC
% Multivariate Granger Causality Toolbox>, _in preparation_, Aug. 2012.
%
%% See also
%
% <autocov_to_smvgc.html |autocov_to_smvgc|> |
% <autocov_to_mvgc_pwc.html |autocov_to_mvgc_pwc|> |
% <autocov_to_var.html |autocov_to_var|> |
% <var2trfun.html |var2trfun|> |
% <autocov_xform.html |autocov_xform|> |
% <sfreqs.html |sfreqs|>
%
%% Copyright notice
%
% _(C) Lionel Barnett, 2012. See file <matlab:edit('license.txt') license.txt> in root directory
% for licensing terms._
%
%%

function [f,fres] = autocov_to_smvgc_pwc_old(G,fres)

[n,~,q1] = size(G);
if nargin < 2 || isempty(fres);
    fres = q1;
end

h = fres+1;
f = nan(n,n,h);

for j = 1:n
    jo = [1:j-1 j+1:n]; % omit j

    [Gj,~,SIGj] = autocov_xform(G,jo); % transform autocov for reduced regression
    [Ajj,SIGjj] = autocov_to_var(Gj);  % transformed VAR parameters
    Hjj = var2trfun(Ajj,fres);         % transformed transfer function

    for ii=1:n-1;
        i  = jo(ii);
        io = [1:ii-1 ii+1:n]; % omit i

        SIGji = SIGjj(io,io)-(SIGjj(io,ii)*SIGjj(ii,io))/SIGjj(ii,ii); % partial covariance
        Hji = Hjj(ii,io,:);                                            % transfer function
        Sji = SIGj(ii,ii);                                             % since reduced residuals serially uncorrelated
        LSji = log(Sji);                                               % spectral density = covariance (flat spectrum)

        for k = 1:h
            f(i,j,k) = LSji - log(real(Sji-Hji(:,:,k)*SIGji*Hji(:,:,k)'));
        end
    end
end
