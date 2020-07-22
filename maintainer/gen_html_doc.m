function gen_html_doc(mname,showCode,stripContents,disHelp)

% Generate html documentation (publish) from source file; mname must be an
% m-file on search path.

if nargin < 2 || isempty(showCode), showCode      = false; end
if nargin < 3,                      stripContents = [];    end
if nargin < 4 || isempty(disHelp),  disHelp       = true;  end

global mvgc_root;

assert(exist(mname,'file') == 2,'can''t find m-file ''%s'' on search path',mname);

popt.showCode = showCode;
popt.evalCode = false;
popt.outputDir = fullfile(mvgc_root,'docs','html');

htmlfile = publish(mname,popt);

fmt_html_doc(mname,stripContents);

fprintf('Generated html file ''%s''\n',htmlfile);

if disHelp
    helpon(mname)
end
