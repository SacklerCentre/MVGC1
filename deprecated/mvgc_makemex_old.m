function [status,result] = mvgc_makemex_old(target)

% Build MVGC 'mex' files, using the makefile in the C folder.

if nargin < 1, target = ''; end

global mvgc_root;

thisdir = cd(fullfile(mvgc_root,'C'));

[status, result] = system(['make MEXEXT=' mexext ' ' target],'-echo');

cd(thisdir);
