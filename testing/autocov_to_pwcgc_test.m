%% autocov_to_pwcgc
%
% Calculate pairwise-conditional time-domain MVGCs (multivariate Granger causalities)
%
% <matlab:open('autocov_to_pwcgc.m') code>
%
%% Syntax
%
%     F = autocov_to_pwcgc(G)
%
%% Arguments
%
% See also <mvgchelp.html#4 Common variable names and data structures>.
%
% _input_
%
%     G          autocovariance sequence
%
% _output_
%
%     F          Granger causality matrix
%
%% Description
%
% Returns the  matrix |F| of pairwise-conditional time-domain MVGCs
%
% <<eq_mvgc_pwc.png>>
%
% (where |[ij]| denotes omission of the |ij|-th variables) between all
% pairs of variables [[ii_inej.png]] represented in |G|, for a stationary VAR
% process with autocovariance sequence |G|. Note that the first index |i| of
% |F| is the target (causee) variable, the second |j| the source (causal)
% variable. See ref. [1] for details.
%
% The caller should take note of any warnings issued by this function and test
% results with a call <isbad.html |isbad|>|(F,false)|.
%
%% References
%
% [1] L. Barnett and A. K. Seth, <matlab:open('mvgc_doc.pdf') The MVGC
% Multivariate Granger Causality Toolbox>, _in preparation_, Aug. 2012.
%
%% See also
%
% <autocov_to_mvgc.html |autocov_to_mvgc|> |
% <autocov_to_var.html |autocov_to_var|> |
% <mvgc_demo.html |mvgc_demo|> |
% <isbad.html |isbad|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function F = autocov_to_pwcgc_test(G)

n = size(G,1);
F = nan(n);

dvtol = 1e-8;

% full regression

owstate = warn_supp;
[~,SIG,DV] = autocov_to_var_test(G);
warn_test(owstate,      'in full regression - bad autocovariance matrix? Check output of ''var_info''');
warn_if(DV(end) > dvtol,'in full regression - regression failed to converge (rel. error = %g)',DV(end));
if warn_if(isbad(SIG),  'in full regression - regression failed'), return; end % show-stopper!

LSIG = log(diag(SIG));

for j = 1:n

    % reduced regression
    
    jo = [1:j-1 j+1:n]; % omit j

    owstate = warn_supp;
    [~,SIGj,DVj] = autocov_to_var_test(G(jo,jo,:));
    warn_test(owstate,       sprintf('in reduced regression for target node %d - bad autocovariance matrix? Check output of ''var_info''',j));
    warn_if(DVj(end) > dvtol,sprintf('in reduced regression for target node %d - regression failed to converge (rel.error = %g)',j,DVj(end)));
    if warn_if(isbad(SIGj),  sprintf('in reduced regression for target node %d - regression failed',j)), return; end % show-stopper!
    
    LSIGj = log(diag(SIGj));

    for ii=1:n-1;
        i = jo(ii);
        F(i,j) = LSIGj(ii)-LSIG(i);
    end
end
