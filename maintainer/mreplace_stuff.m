function mreplace_stuff(mdir)

if nargin < 1 || isempty(mdir)
    mreplace_stuff('core');
    mreplace_stuff(['core' filesep 'GCCA_compat']);
    mreplace_stuff(['core' filesep 'subsample']);
    mreplace_stuff('demo');
    mreplace_stuff('docs');
    mreplace_stuff('experimental');
    mreplace_stuff('gc');
    mreplace_stuff(['gc' filesep 'GCCA_compat']);
    mreplace_stuff(['gc' filesep 'subsample']);
    mreplace_stuff('stats');
    mreplace_stuff('utils');
    mreplace_stuff(['utils' filesep 'control']);
    mreplace_stuff(['utils' filesep 'legacy']);
    mreplace_stuff(['utils' filesep 'legacy' filesep 'randi']);
    mreplace_stuff(['utils' filesep 'legacy' filesep 'rng']);
    mreplace_stuff(['utils' filesep 'stats']);
    
    % maybe not these always
    % mreplace_stuff(['docs' filesep 'html']);
    % mreplace_stuff('deprecated');
    % mreplace_stuff('testing');
    return
end

% Replace stuff in m files

global mvgc_root;

mdir = fullfile(mvgc_root,mdir);

mcells = struct2cell(dir([mdir filesep '*.m']));
n = size(mcells,2); % 1st two are ',' and '..'
for i = 1:n
    mfile = mcells{1,i};
    fprintf('Processing m file %2d of %d = ''%s'' : ',i,n,mfile);
    
    mfilef = fullfile(mdir,mfile);
    
    % back it up!
    
    syscmd = ['cp ' mfilef ' ' tempdir 'mreplace_stuff' filesep];
    status = system(syscmd,'-echo');
    if status == 0
        fprintf('backed up ok : ');
    else
       fprintf(2,'failed to backup file ''%s''\n',mfilef);
       continue
    end

    fidr = fopen(mfilef, 'r');
    if fidr == -1
        fprintf(2,'failed to open input file ''%s''\n',mfilef);
        continue
    end    
    str = fread(fidr,inf,'*char')';    
    status = fclose(fidr);
    if status ~=  0
        fprintf(2,'failed to close input file ''%s''\n',mfilef);
        continue
    end
    
    %%%   DO THE BUSINESS    %%%
%{
    expr   = 'matlab:edit(';
    repstr = 'matlab:open(';
    str    = strrep(str,expr,repstr);
    
    expr   = 'matlab:openpdf(';
    repstr = 'matlab:open(';
    str    = strrep(str,expr,repstr);
    
    expr   = sprintf('%% _(C) Lionel Barnett, 2012. See file <matlab:open(''license.txt'')\n%% license.txt> in root directory for licensing terms._');
    repstr = sprintf('%% [(C)] _Lionel Barnett and Anil K. Seth, 2012. See file\n%% <matlab:open(''license.txt'') license.txt> in root directory for licensing\n%% terms._');
    str    = strrep(str,expr,repstr);
    
    expr   = sprintf('%% [(C)] _Lionel Barnett and Anil K. Seth, 2012. See file\n%% <matlab:open(''license.txt'') license.txt> in root directory for licensing\n%% terms._');
    repstr = sprintf('%% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in\n%% installation directory for licensing terms.');
    str    = strrep(str,expr,repstr);
    
    expr   = '%% Copyright(.*?)% \(C\)';
    repstr = '% \(C\)';
    str    = regexprep(str,expr,repstr);
    
    expr   = 'autocov_to_mvgc_pwc';
    repstr = 'autocov_to_pwcgc';
    str    = strrep(str,expr,repstr);
    
    expr   = 'autocov_to_smvgc_pwc';
    repstr = 'autocov_to_spwcgc';
    str    = strrep(str,expr,repstr);
    
    expr   = 'GCCA_tsdata_to_mvgc_pwc';
    repstr = 'GCCA_tsdata_to_pwcgc';
    str    = strrep(str,expr,repstr);
    
    expr   = '<matlab:open(''mvgc_doc.pdf'')';
    repstr = '<matlab:open(''mvgc.pdf'')';
    str    = strrep(str,expr,repstr);
    
    expr   = 'Multivariate Granger Causality Toolbox>, _in preparation_, Aug. 2012.';
    repstr = sprintf('Multivariate Granger Causality Toolbox: A New Approach to Granger-causal\n%% Inference>, _J. Neurosci. Methods_ 223, 2014.');
    str    = strrep(str,expr,repstr);
    
    expr   = '_J. Neurosci. Methods_, 223, 2014';
    repstr = '_J. Neurosci. Methods_ 223, 2014';
    str    = strrep(str,expr,repstr);
%}
    
    expr   = sprintf('%% [1] L. Barnett and A. K. Seth, <matlab:open(''mvgc_doc.pdf'') The MVGC\n%% Multivariate Granger Causality Toolbox: A New Approach to Granger-causal\n%% Inference>, _J. Neurosci. Methods_ 223, 2014.');
    repstr = sprintf('%% [1] L. Barnett and A. K. Seth,\n%% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC\n%%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal\n%% Inference>, _J. Neurosci. Methods_ 223, 2014\n%% [ <matlab:open(''mvgc_preprint.pdf'') preprint> ].');
    str    = strrep(str,expr,repstr);
    
    %%% THE BUSINESS IS DONE %%%

    fidw = fopen(mfilef, 'w');
    if fidw == -1
        fprintf(2,'failed to open output file ''%s''\n',mfilef);
        continue
    end
    fwrite(fidw,str);
    status = fclose(fidw);
    if status ~=  0
        fprintf(2,'failed to close output file ''%s''\n',mfilef);
        continue
    end
 
    fprintf('done\n');
end
