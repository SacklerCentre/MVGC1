%% plot_autocov
%
% Autocovariance plotting utility
%
% <matlab:open('plot_autocov.m') code>
%
%% Syntax
%
%     plot_autocov(G,cflag)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     G          autocovariance sequence
%     cflag      plot autocorrelation rather than autocovariance (default: 0)
%
%% Description
%
% Plots autocovariance sequence |G| against lags for each pair of variables on a
% grid. If the |cflag| flag is nonzero, autocorrelation rather than
% autocovariance is plotted (see <cov2corr.html |cov2corr|>). If |cflag| is
% negative, then zero-lag correlation is omitted (this may be useful since
% lag-zero autocorrelation can be an order of magnitude larger than
% autocorrelation at higher lags - e.g. for a near-white noise process).
%
%% See also
%
% <cov2corr.html |cov2corr|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function plot_autocov_old(G,cflag)

if nargin < 2 || isempty(cflag), cflag = 0; end

[n,n1,q1] = size(G);
assert(n1 == n,'autocovariance matrix has bad shape');
q = q1-1;

if cflag ~= 0
    G = cov2corr(G);
    if cflag > 0
        t = (0:q)';
        xlims = [0 q];
        ylims = [-1 1];
    else
        G = G(:,:,2:end); % omit zero-lag correlation
        t = (1:q)';
        xlims = [1 q];
        ylims = [min(G(:)) max(G(:))];
    end
else
    t = (0:q)';
    xlims = [0 q];
    ylims = [min(G(:)) max(G(:))];
end

k = 0;
for i = 1:n
    for j = 1:n
        k = k+1;
        subplot(n,n,k);
        plot(t,squeeze(G(i,j,:)));
        grid on;
        xlabel('lags');
        ylabel(sprintf('var %d,%d',i,j));
        xlim(xlims);
        ylim(ylims);
    end
end
