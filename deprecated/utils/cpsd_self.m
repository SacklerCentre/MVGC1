function SD = cpsd_self(S)

% SD = cpsd_self(S)
%
% Return matrix of self-power spectra (without cross-spectra).

[n,~,h] = size(S);

SD = zeros(n,h);
for k = 1:h
    SD(:,k) = diag(S(:,:,k));
end
