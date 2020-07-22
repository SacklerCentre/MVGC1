function regen_html_docs(mdir,showCode,stripContents)

% Regenerate html documentation (publish) from source files in directory
% mdir

global mvgc_root;

if nargin < 1 || isempty(mdir)
    regen_html_docs('.',true,false);
    regen_html_docs('core',false,true);
    regen_html_docs('gc',false,true);
    regen_html_docs(['gc' filesep 'GCCA_compat'],false,true);
    regen_html_docs(['gc' filesep 'subsample'],false,true);
    regen_html_docs('demo',true,false);
    regen_html_docs('docs',true,false);
    regen_html_docs('experimental',false,true);
    regen_html_docs('stats',false,true);
    regen_html_docs('utils',false,true);
    
    if ~fexists(@rng) || ~fexists(@randi) % legacy hack
        regen_html_docs(['utils' filesep 'legacy'],false,true);
        if ~fexists(@randi)
            regen_html_docs(['utils' filesep 'legacy' filesep 'randi'],false,true);
        else
            addpath(fullfile(mvgc_root,'utils','legacy','randi'));
            regen_html_docs(['utils' filesep 'legacy' filesep 'randi'],false,true);
            rmpath(fullfile(mvgc_root,'utils','legacy','randi'));
        end
        if ~fexists(@rng)
            regen_html_docs(['utils' filesep 'legacy' filesep 'rng'],false,true);
        else
            addpath(fullfile(mvgc_root,'utils','legacy','rng'));
            regen_html_docs(['utils' filesep 'legacy' filesep 'rng'],false,true);
            rmpath(fullfile(mvgc_root,'utils','legacy','rng'));
        end
    else
        addpath(fullfile(mvgc_root,'utils','legacy'));
        regen_html_docs(['utils' filesep 'legacy'],false,true);
        rmpath(fullfile(mvgc_root,'utils','legacy'));
    end
    if ~fexists(@dlyap) % Control System Toolbox hack
        regen_html_docs(['utils' filesep 'control'],false,true);
    else
        addpath(fullfile(mvgc_root,'utils','control'));
        regen_html_docs(['utils' filesep 'control'],false,true);
        rmpath(fullfile(mvgc_root,'utils','control'));
    end
    if ~fexists(@chi2cdf) % Statistics Toolbox hack
        regen_html_docs(['utils' filesep 'stats'],false,true);
    else
        addpath(fullfile(mvgc_root,'utils','stats'));
        regen_html_docs(['utils' filesep 'stats'],false,true);
        rmpath(fullfile(mvgc_root,'utils','stats'));
    end
    
    return
end

if nargin < 2 || isempty(showCode), showCode      = false; end
if nargin < 3,                      stripContents = [];    end

mdir = fullfile(mvgc_root,mdir);

mcells = struct2cell(dir(mdir));
n = size(mcells,2)-2; % 1st two are ',' and '..'
for i = 1:n
    [pathstr,name,ext] = fileparts(mcells{1,i+2});
    if strcmp(ext,'.m')
        mfile = fullfile(pathstr,name); % strip extension
        gen_html_doc(mfile,showCode,stripContents,false);
    end
end
