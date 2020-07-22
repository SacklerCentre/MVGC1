n         = 4;      % number of variables
p         = 3;      % autoregressive order

rho       = 0.9;    % AR spectral norm
icfac     = 10;     % noise cross-correlation factor (bigger = less correlated, Inf for completely uncorrelated)
mseed     = 0;      % model random seed (0 for unseeded)

fres      = 1000;   % frequency resolution

numtol    = 1e-8;   % numerical tolerance for convergence

%-------------------------------------------------------------------------------

% generate random VAR parameters

rstate = rng_seed(mseed);
A = var_specrad(randn(n,n,p),rho);
V = cov_rand(n,icfac);
rng_restore(rstate);

% calculate exact CPSD and transfer function for VAR. (NOTE: we use a new
% var2cpsd routine rather than the MVGC v1.0 var_to_cpsd, since the latter
% appears to produce slighlty off-Hermitian S, which causes wilson_sf to error
% out.

[S,H] = var2cpsd(A,V,fres);

% MVGC Wilson implementation

fprintf('cpsd_to_var : ');
tic
[H1,V1] = cpsd_to_var(S,[],[],numtol);
S1 = HV2S(H1,V1);
toc

% wilson_sf

fprintf('wilson_sf   : ');
tic
[H2,V2] = wilson_sf(S,1,numtol);
S2 = HV2S(H2,V2);
toc

% report relative errors for CPSD, transfer function and noise covariance matrix

Serr1 = norm(S1(:)-S(:),2)/norm(S(:),2);
Herr1 = norm(H1(:)-H(:),2)/norm(H(:),2);
Verr1 = norm(V1(:)-V(:),2)/norm(V(:),2);
fprintf('\nrel. errors (cpsd_to_var) : S = %.2e, H = %.2e, V = %.2e\n',Serr1,Herr1,Verr1);

Serr2 = norm(S2(:)-S(:),2)/norm(S(:),2);
Herr2 = norm(H2(:)-H(:),2)/norm(H(:),2);
Verr2 = norm(V2(:)-V(:),2)/norm(V(:),2);
fprintf('rel. errors (wilson_sf)   : S = %.2e, H = %.2e, V = %.2e\n',Serr2,Herr2,Verr2);
