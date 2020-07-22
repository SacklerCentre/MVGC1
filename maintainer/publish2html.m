function publish2html(www_dir,htmlfiles,copyStuff,procHtml)

% Publish html documentation

global mvgc_root;
global mvgc_version;

if nargin < 1 || isempty(www_dir),   www_dir   = fullfile(tempdir,'MVGC'); end
if nargin < 2,                       htmlfiles = [];                       end
if nargin < 3 || isempty(copyStuff), copyStuff = true;                     end
if nargin < 3 || isempty(procHtml),  procHtml  = true;                     end

verstr = ['v' num2str(mvgc_version.major) '.' num2str(mvgc_version.minor)];

www_html_dir = fullfile(www_dir,'html');
www_imag_dir = fullfile(www_dir,'images');
www_docs_dir = fullfile(www_dir,'docs');

makedir(www_dir);
makedir(www_html_dir);
makedir(www_docs_dir);
makedir(www_imag_dir);

src_dir            = mvgc_root;
src_maintainer_dir = fullfile(src_dir,'maintainer');
src_docs_dir       = fullfile(src_dir,'docs');
src_docs_html_dir  = fullfile(src_docs_dir,'html');

if copyStuff
    copyfiles(src_maintainer_dir, www_dir,      'index.html');
    copyfiles(src_dir,            www_dir,      'license.txt');
    copyfiles(src_docs_dir,       www_docs_dir, 'mvgc_preprint.pdf');
    copyfiles(src_docs_html_dir,  www_imag_dir, '*.png');
end

if ~procHtml, return; end

% process all html docs

if isempty(htmlfiles)
    htmlcells = struct2cell(dir([src_docs_html_dir filesep '*.html']));
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
    fprintf('Processing html file %2d of %d = ''%s'' : ',i,n,htmlfile);
    
    infile  = fullfile(src_docs_html_dir,htmlfile);
    outfile = fullfile(www_html_dir,htmlfile);
    
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
    
    fnames = {'paretotails','kstest'};
    for j = 1:length(fnames)
        expr   = ['matlab/ref/' fnames{j}];
        repstr = ['stats/' fnames{j}];
        str    = strrep(str,expr,repstr);
    end
     
    % license
    
    expr   = '<a href="matlab:open(''license.txt'')">license.txt</a> in installation directory';
    repstr = '<a href="../license.txt">license.txt</a>';
    str    = strrep(str,expr,repstr);
    
    % mvgc_doc
    
    expr   = '<a href="matlab:open(''mvgc_preprint.pdf'')">';
    repstr = '<a href="../docs/mvgc_preprint.pdf">';
    str    = strrep(str,expr,repstr);
    
    % code
    
    expr   = '<p><a href="matlab:open\((.*?)\)">code</a></p>';
    repstr = '';
    str    = regexprep(str,expr,repstr);
    
    % images
    
    expr   = 'src="(.*?)\.png"';
    repstr = 'src="../images/$1.png"';
    str    = regexprep(str,expr,repstr);
    
    % strip source
    
    expr   = '##### SOURCE BEGIN #####(.*?)##### SOURCE END #####';
    repstr = '';
    str    = regexprep(str,expr,repstr);
    
    % munge email address to avoid spambots
    
    expr   = '<a href="mailto:mvgctoolbox@sussex.ac.uk">mvgctoolbox @ sussex.ac.uk</a>';
    repstr = 'm v g c t o o l b o x AT s u s s e x . a c . u k';
    str    = strrep(str,expr,repstr);

    %%% THE BUSINESS IS DONE %%%
    
    fidw = fopen(outfile, 'w');
    if fidw == -1
        fprintf(2,'failed to open output file ''%s''\n',outfile);
        continue
    end

    fwrite(fidw,str);
    
    status = fclose(fidw);
    if status ~=  0
        fprintf(2,'failed to close output file ''%s''\n',outfile);
        continue
    end
    
    fprintf('done\n');
end

function makedir(dirname)

syscmd = ['mkdir -p ' dirname];
%fprintf('\nmake dir command ''%s'':\n\n',syscmd);
status = system(syscmd,'-echo');
if status ~= 0
   fprintf(2,'command ''%s'' failed\n',syscmd);
end

function copyfiles(sdir,tdir,filename)

syscmd = ['cp -v ' sdir filesep filename ' ' tdir filesep];
%fprintf('\ncopy command ''%s'':\n\n',syscmd);
fprintf('\n');
status = system(syscmd,'-echo');
if status ~= 0
   fprintf(2,'command ''%s'' failed\n',syscmd);
end
