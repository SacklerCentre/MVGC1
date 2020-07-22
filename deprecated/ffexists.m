function res = ffexists(fname)

try
    eval(['[~] = ' fname]); % the [~] is to suppress output
    %fprintf('no exception\n');
    res = true;
catch except
    %fprintf('exception: ''%s'', ''%s''\n', except.identifier,except.message);
    res = ~strcmpi(except.identifier,'MATLAB:UndefinedFunction');
end
