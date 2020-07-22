function publish2html_old1(tdir,htmlfiles,copyImages,copyLicense,linkIndex)

% Publish html documentation

global mvgc_root;

if nargin < 1 || isempty(tdir), tdir = fullfile(tempdir,'MVGC','htmldocs'); end

if nargin < 2, htmlfiles = [];  end

if nargin < 3 || isempty(copyImages),  copyImages  = false; end
if nargin < 4 || isempty(copyLicense), copyLicense = false; end
if nargin < 4 || isempty(linkIndex),   linkIndex   = false; end

sdir = fullfile(mvgc_root,'docs','html');

if isempty(htmlfiles)
    htmlcells =  struct2cell(dir([sdir filesep '/*.html']));
    n = size(htmlcells,2);
    htmlfiles = cell(n,1);
    for i = 1:n
        htmlfiles{i} = htmlcells{1,i};
    end
else
    if iscell(htmlfiles)
        n = length(htmlfiles);
    else
        n = 1;
        tmp = cell(1,1);
        tmp{1} = htmlfiles;
        htmlfiles = tmp;
    end        
    for i = 1:n
        htmlfiles{i} = [htmlfiles{i} '.html'];
    end
end

for i = 1:n
    htmlfile = htmlfiles{i};
    fprintf('Processing html file %2d of %d =''%s'' : ',i,n,htmlfile);
    
    infile  = fullfile(sdir,htmlfile);
    outfile = fullfile(tdir,htmlfile);
    
    fidr = fopen(infile, 'r');
    if fidr == -1
        fprintf(2,'failed to open input file ''%s''\n',infile);
        continue
    end
    
    str = fread(fidr,inf,'*char')';
    
    status = fclose(fidr);
    if status ~=  0
        fprintf(2,'failed to close input file ''%s''\n',infile);
        continue
    end
    
    fidw = fopen(outfile, 'w');
    if fidw == -1
        fprintf(2,'failed to open output file ''%s''\n',outfile);
        continue
    end
    
    %%%   DO THE BUSINESS    %%%
    
    % matlab:doc
    
    expr   = '<a href="matlab:doc\(''(.*?)''\)">';
    repstr = '<a href="http://www.mathworks.com/help/matlab/ref/$1.html">';
    str    = regexprep(str,expr,repstr);
    
    expr   = 'matlab/ref/dlyap';
    repstr = 'control/ref/dlyap';
    str    = strrep(str,expr,repstr);
    
    fnames = {'pwelch','cpsd','dpss'};
    for j = 1:length(fnames)
        expr   = ['matlab/ref/' fnames{j}];
        repstr = ['signal/ref/' fnames{j}];
        str    = strrep(str,expr,repstr);
    end
    
    expr   = 'matlab/ref/paretotails';
    repstr = 'stats/paretotails';
    str    = strrep(str,expr,repstr);
    
    % license
    
    expr   = '<a href="matlab:open(''license.txt'')">license.txt</a> in root directory';
    repstr = '<a href="license.txt">license.txt</a>';
    str    = strrep(str,expr,repstr);
    
    % code
    
    expr   = '<p><a href="matlab:open\((.*?)\)">code</a></p>';
    repstr = '';
    str    = regexprep(str,expr,repstr);
    
    %%% THE BUSINESS IS DONE %%%

    fwrite(fidw,str);
    
    status = fclose(fidw);
    if status ~=  0
        fprintf(2,'failed to close output file ''%s''\n',outfile);
        continue
    end
    
    fprintf('done\n');
end

% copy images across

if copyImages
    syscmd = ['cp -v ' sdir filesep '*.png ' tdir filesep];
    fprintf('\nGoing for image copy with command ''%s'':\n\n',syscmd);
    status = system(syscmd,'-echo');
    if status ~= 0
       fprintf(2,'command ''%s'' failed\n',syscmd);
    end
end

% copy license file across

if copyLicense
    syscmd = ['cp -v ' mvgc_root filesep 'license.txt ' tdir filesep];
    fprintf('\nGoing for license copy with command ''%s'':\n\n',syscmd);
    status = system(syscmd,'-echo');
    if status ~= 0
       fprintf(2,'command ''%s'' failed\n',syscmd);
    end
end

% link index

if linkIndex
    syscmd = ['ln -s ' tdir filesep 'mvgcindex.html ' tdir filesep 'index.html'];
    fprintf('\nGoing for index link with command ''%s'':\n\n',syscmd);
    status = system(syscmd,'-echo');
    if status ~= 0
       fprintf(2,'command ''%s'' failed\n',syscmd);
    end
end
