function [AF,VF,AA] = autocov_to_var_diag(G)

tol = 1e-16;

[n,~,q1] = size(G);
q = q1-1;
qn = q*n;

G0 = G(:,:,1);                                               % covariance
GF = reshape(G(:,:,2:end),n,qn)';                            % forward  autocov sequence
GB = reshape(permute(flipdim(G(:,:,2:end),3),[1 3 2]),qn,n); % backward autocov sequence

AF = zeros(n,qn); % forward  coefficients
AB = zeros(n,qn); % backward coefficients (reversed compared with Whittle's treatment)
VF = zeros(n,n,q);

AA = zeros(n,n,q);

% initialise recursion

k = 1;            % model order

r = q-k;
kf = 1:k*n;       % forward  indices
kb = r*n+1:qn;    % backward indices

AF(:,kf) = GB(kb,:)/G0;
AB(:,kb) = GF(kf,:)/G0;

% and loop

VF(:,:,1) = G0-AF(:,kf)*GF(kf,:);

for k=2:q

    AAF = (GB((r-1)*n+1:r*n,:)-AF(:,kf)*GB(kb,:))/(G0-AB(:,kb)*GB(kb,:)); % DF/VB
    AAB = (GF((k-1)*n+1:k*n,:)-AB(:,kb)*GF(kf,:))/(G0-AF(:,kf)*GF(kf,:)); % DB/VF

    AFPREV = AF(:,kf);
    ABPREV = AB(:,kb);

    r = q-k;
    kf = 1:k*n;
    kb = r*n+1:qn;

    AF(:,kf) = [AFPREV-AAF*ABPREV AAF];
    AB(:,kb) = [AAB ABPREV-AAB*AFPREV];
    
    VF(:,:,k) = G0-AF(:,kf)*GF(kf,:);
    AA(:,:,k) = reshape(AF(:,(k-1)*n+1:k*n),n,n);
    % rdiff = abs(det(VF(:,:,k))-det(VF(:,:,k-1)))/det(VF(:,:,k-1));
    % rdiff = abs(norm(VF(:,:,k)-VF(:,:,k-1)))/norm(VF(:,:,k-1));
    % if rdiff < tol, break; end
end

VF = VF(:,:,1:k);
AF = reshape(AF(:,1:k*n),n,n,k);
