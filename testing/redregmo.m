S      = 1000;

m      = 10000;      % number of observations

tstat     = '';     % statistical test for MVGC:  'chi2' for Geweke's chi2 test (default) or'F' for Granger's F-test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

fs        = 200;    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)

seed      = 123872; % random seed (0 for unseeded)

rho  =  0.98;

a    =  0.6;
c    = -0.7;
b    =  0.5;

mof  = 1;
mor  = 2;

nocause = 1;

%-------------------------------------------------------------------------------

rng_seed(seed);

n = 2;
x = 1;
y = 2;

A = [a b; 0 c];
A = var_specrad(A,rho);

if nocause, A(1,2) = 0; end

a = A(1,1);
b = A(1,2);
c = A(2,2);

V = eye(n);

[G,res] = var_to_autocov(A,V);
fprintf('\nVAR check:\n'); disp(res);
assert(~res.error,'bad VAR');

%FTi = autocov_to_mvgc(G,x,y);
DELTA = (1+b*b+c*c)/2;
FTi   = log(DELTA + sqrt(DELTA*DELTA-c*c));

aclags = size(G,3)-1;
Gp = G(:,:,1:min(aclags,mor)+1);
FTp = autocov_to_mvgc(Gp,x,y);
%   p == 1 case
%{
gamma0 = 1/(1-c*c);
beta0  = ((b*c)/(1-a*c))*gamma0;
alpha0 = (1+2*a*b*beta0+b*b*gamma0)/(1-a*a);
alpha1 = a*alpha0+b*beta0;
FTp    = log(alpha0-(alpha1*alpha1)/alpha0);
%}

Fi = zeros(S,1);
Fp = zeros(S,1);
Fq = zeros(S,1);

s = 0;
while s < S;
    fprintf('sample %d of %d',s+1,S);
    
    U = var_to_tsdata(A,V,m);

    [Ae,Ve] = tsdata_to_var(U,mof,'OLS');

    [Ge,rese] = var_to_autocov(Ae,Ve);
    if rese.error
        fprintf(2,' - bad VAR: %s\n',rese.errmsg);
        continue;
    end

    s = s+1;
    
    Fi(s) = autocov_to_mvgc(Ge,x,y);
    
    Fp(s) = GCCA_tsdata_to_mvgc(U,x,y,[mof,mor],'OLS');
    
    Fq(s) = GCCA_tsdata_to_mvgc(U,x,y,mor,'OLS');
    
    fprintf('\n');
end

fprintf('\nFTi = %.8f\n',FTi);
fprintf('Fi  = %.8f\n',mean(Fi));
fprintf('\nFT1 = %.8f\n',FTp);
fprintf('Fp  = %.8f\n',mean(Fp));
fprintf('Fq  = %.8f\n',mean(Fq));
