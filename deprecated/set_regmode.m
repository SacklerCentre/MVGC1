function oldmode = set_regmode(rm)

% Set regression mode

global REGMODE;

oldmode = REGMODE;

if nargin < 1 || isempty(rm), rm = 1; end % default is OLS

errstr = 'regression mode must be a string (''OLS'' or ''LWR'') or a number (1 or 2)';

if isscalar(rm) && isnumeric(rm)
    assert(rm >= 1 && rm <= 2,errstr);
    REGMODE = rm;
else
    assert(ischar(rm),errstr);
    switch upper(rm)
        case 'OLS', REGMODE = 1;
        case 'LWR', REGMODE = 2;
        otherwise, error(errstr);
    end
end

if nargout == 0
    get_regmode;
end
