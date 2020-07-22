function [rmode,rmstr] = get_regmode

% Get regression mode

global REGMODE;

switch REGMODE
    case 1, rmstr = 'OLS';
    case 2, rmstr = 'LWR';
    otherwise, error('unknown regression mode!');
end

if nargout == 0
    fprintf('regression mode is set to %d (%s)\n',REGMODE,rmstr);
else
    rmode = REGMODE;
end
