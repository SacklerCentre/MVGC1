function [status, result] = openpdf(filename)

% Open pdf document (roll your own)

assert(exist(filename, 'file') == 2,'could not find file ''%s''',filename);
fname = which(filename);
if isempty(fname)
    fname = filename;
end
fprintf('opening file ''%s''\n', fname);
[status, result] = system(['xpdf ' fname ' &'],'-echo');
