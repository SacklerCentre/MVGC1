function S = HV2S(H,V)

[n,~,h] = size(H);
S = NaN;
[VSR,cholp] = chol(V,'lower');
if cholp, return; end % show stopper
S = zeros(n,n,h);
for k = 1:h
    HVSRk = H(:,:,k)*VSR;
    S(:,:,k) = HVSRk*HVSRk';
end
