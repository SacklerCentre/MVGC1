function AF = var_scoef(A,fres)

% AF = var_scoef(A,fres)
%
% Return spectral VAR coefficients AF for VAR with coefficients A. fres
% specifies the frequency resolution; call 'sfreqs(fres,fs)' where fs is
% the sample rate to get a corresponding vector of frequencies on [0,fs/2].

n = size(A,1);
I = eye(n);
AF = bfft(cat(3,I,-A),2*fres); % over [0,2*pi)
