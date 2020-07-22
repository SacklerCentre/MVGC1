%% plot_cpsd
%
% Cross-power spectral density plotting utility
%
% <matlab:open('plot_cpsd.m') code>
%
%% Syntax
%
%     plot_cpsd(lam,auto,S1,lstr1,S2,lstr2,...)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     lam        vector of frequencies
%     auto       auto-spectra only (default: true)
%     Sn         n-th cpsd matrix
%     lstrn      legend string for n-th cpsd plot
%
%% Description
%
% Plots multiple cpsds vs. frequency on a grid. (|Sn|, |lstrn|) specify
% (cpsd, legend) pairs. If the |auto| flag is set only auto-spectra are
% plotted.
%
%% See also
%
% <mvgc_demo_stats.html |mvgc_demo_stats|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function plot_cpsd_old(auto,lam,varargin)

if isempty(auto), auto = true; end % default is auto-spectra only

assert(isvector(lam),'frequencies must be a vector');
h = length(lam);
lam = lam(:);

N = length(varargin)/2;
assert(2*N == length(varargin) && N > 0,'bad number of arguments');
[n,n1,h1] = size(varargin{1});
assert(n1 == n,'cpsd matrix 1 has bad shape');
assert(h1 == h,'cpsd matrix does not match frequency vector');

S = complex(zeros(n,n,h,N));
leg = cell(N,1);
for s = 1:N
    assert(isequal([n n h],size(varargin{2*s-1})),'cpsd matrix %d does not match cpsd matrix 1',s);
    assert(ischar(varargin{2*s}),'legend %d not a string',s);
    S(:,:,:,s) = varargin{2*s-1};
    leg{s} = varargin{2*s};
end

if auto % auto-spectra only

    for i = 1:n
        subplot(n,1,i);
        plot(lam,squeeze(S(i,i,:,:)));
        xlim([lam(1),lam(end)]);
        xlabel('frequency');
        ylabel(sprintf('var %d auto-power',i));
        legend(leg);
    end

else    % auto- and cross-spectra

    k = 0;
    for i = 1:n
        for j = 1:n
            k = k+1;
            if i == j
                subplot(n,n,k);
                plot(lam,squeeze(S(i,i,:,:)));
                xlim([lam(1),lam(end)]);
                xlabel('frequency');
                ylabel(sprintf('var %d auto-power',i));
                legend(leg);
            elseif i < j
                subplot(n,n,k);
                plot(lam,abs(squeeze(S(i,j,:,:))));
                xlim([lam(1),lam(end)]);
                xlabel('frequency');
                ylabel(sprintf('var %d,%d cross-power',i,j));
                legend(leg);
            end
        end
    end

end
