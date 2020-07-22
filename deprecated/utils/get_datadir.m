function datadir = get_datadir(subdir)

datadir = getenv('DATADIR');
assert(~isempty(datadir),'environment variable ''DATADIR'' not found');
if nargin > 0 && ~isempty(subdir)
    datadir = fullfile(datadir,subdir);
end
if exist(datadir,'dir') ~= 7
    fprintf(2,'WARNING (get_datadir): data directory ''%s'' not found\n',datadir);
end
