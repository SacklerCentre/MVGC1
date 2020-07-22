function [status, result] = openpng(name)

% Open pdf document (roll your own)

assert(exist(name,'file') == 2,'can''t find ''%s'' on search path',name);
fullname = which(name);
[status, result] = system(['xdg-open ' fullname ' &'],'-echo');
