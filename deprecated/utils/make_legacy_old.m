function make_legacy_old(tdir)

global mvgc_root;

if nargin < 1 || isempty(tdir), tdir = [tempdir 'mvgc_legacy']; end

fprintf('Populating target directory...');
syscmd = ['cp -r ' mvgc_root ' ' tdir];
status = system(syscmd,'-echo');
if status == 0
    fprintf(' done\n');
else
   fprintf(2,' failed\n');
   return
end

make_leg([],tdir);
make_leg('core',tdir);
make_leg(['core' filesep 'GCCA_compat'],tdir);
make_leg(['core' filesep 'subsample'],tdir);
make_leg('demo',tdir);
make_leg('docs',tdir);
make_leg('experimental',tdir);
make_leg('maintainer',tdir);
make_leg('stats',tdir);
make_leg('testing',tdir);
make_leg('utils',tdir);
make_leg(['utils' filesep 'legacy'],tdir);
make_leg(['utils' filesep 'legacy' filesep 'randi'],tdir);
make_leg(['utils' filesep 'legacy' filesep 'rng'],tdir);

function make_leg(mdir1,tdir)

global mvgc_root;

% Replace stuff in m files

mdir = fullfile(mvgc_root,mdir1);
tdir = fullfile(tdir,mdir1);

fprintf('Target ''%s''\n',tdir);

mcells = struct2cell(dir([mdir filesep '*.m']));
n = size(mcells,2); % 1st two are ',' and '..'
for i = 1:n
    mfile = mcells{1,i};
    fprintf('\tProcessing m file %2d of %d = ''%s''...',i,n,mfile);
    
    mfilef = fullfile(mdir,mfile);
    fidr = fopen(mfilef, 'r');
    if fidr == -1
        fprintf(2,' failed to open input file ''%s''\n',mfilef);
        continue
    end    
    str = fread(fidr,inf,'*char')';    
    status = fclose(fidr);
    if status ~=  0
        fprintf(2,' failed to close input file ''%s''\n',mfilef);
        continue
    end
    
    % replace '~' in return values with 'dummy'

    expr   = '[~,';
    repstr = '[dummy,';
    str    = strrep(str,expr,repstr);

    expr   = ',~,';
    repstr = ',dummy,';
    str    = strrep(str,expr,repstr);

    expr   = ',~]';
    repstr = ',dummy]';
    str    = strrep(str,expr,repstr);

    expr   = '[~]';
    repstr = 'dummy';
    str    = strrep(str,expr,repstr);

    tfilef = fullfile(tdir,mfile);
    fidw = fopen(tfilef, 'w');
    if fidw == -1
        fprintf(2,' failed to open output file ''%s''\n',tfilef);
        continue
    end
    fwrite(fidw,str);
    status = fclose(fidw);
    if status ~=  0
        fprintf(2,' failed to close output file ''%s''\n',tfilef);
        continue
    end
 
    fprintf(' done\n');
end
