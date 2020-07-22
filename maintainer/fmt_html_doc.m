function str = fmt_html_doc(mname,stripContents)

% Format (published) html documentation

if nargin < 2 || isempty(stripContents), stripContents = true;  end

global mvgc_root;
global mvgc_version;

verstr = ['v' num2str(mvgc_version.major) '.' num2str(mvgc_version.minor)];

htmlfile = fullfile(mvgc_root,'docs','html',[mname '.html']);

assert(exist(htmlfile,'file') == 2,'can''t find html help file ''%s''',htmlfile);

fid = fopen(htmlfile, 'r');
assert(fid ~= -1,'failed to open html help file ''%s'' for reading',htmlfile);
str = fread(fid,inf,'*char')';
status = fclose(fid);
assert(status ==  0,'failed to close html help file ''%s''',htmlfile);

if stripContents
    str = regexprep(str, '<h2>Contents</h2><div>(.*?)</div>', '');
end

% indent equations

str = strrep(str, 'img vspace="5" hspace="5" src="eq_', 'img vspace="5" hspace="24" src="eq_');

% inline images

str = regexprep(str, '\[\[(.*?)\]\]', '<img valign="middle" src="$1">');

% other stuff

str = strrep(str, '[R^2]', 'R<sup>2</sup>');
str = strrep(str, '[chi^2]', '&chi;<sup>2</sup>');
str = strrep(str, '[lambda]', '&lambda;');
str = strrep(str, '[H_0]', 'H<sub>0</sub>');
%str = strrep(str, '[(C)]', '&copy;');

expr = '<p>\(C\) Lionel Barnett(.*?)</p>';
repstr = '';
str = regexprep(str,expr,repstr);

mver = version('-release');
if strcmpi(mver,'2011a')
    expr = 'Published with MATLAB&reg; 7.12';
else
    expr = ['<a href="http://www.mathworks.com/products/matlab/">Published with MATLAB&reg; R' version('-release')];
end
repstr = ['MVGC Toolbox ' verstr '. &copy; Lionel Barnett and Anil K. Seth, 2012.<br>See file <a href="matlab:open(''license.txt'')">license.txt</a> in installation directory for licensing terms.'];
str = strrep(str,expr,repstr);


fid = fopen(htmlfile, 'w');
assert(fid ~= -1,'failed to open find html help file ''%s'' for writing',htmlfile);
fwrite(fid,str,'char');
status = fclose(fid);
assert(status ==  0,'failed to close find html help file ''%s''',htmlfile);
