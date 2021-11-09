%% sgtitlex
%
% Call |sgtitle| if available, else do nothing
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function sgtitlex(ptitle)

if exist('sgtitle');
	sgtitle(ptitle);
% else
%	do nothing
end
