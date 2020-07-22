function [AX,XZ,EE] = freedman_lane(U,p,x,y)

% [AX,XZ,EE] = freeman_lane(U,p,x,y)
%
% Return critical variables for generation of permutation data according to
% Freedman and Lane method for a null hypothesis of no dependence (up to p
% lags) of X on Y, conditional on remaining variables.
%
% To generate a permutation UP, the sequence is
%
%    [AX,XZ,EE] = freedman_lane(U,p,x,y)
%    UP = U(:,p+1:m);
%    UP(x,:) = genvar(AX,XZ+EE(:,randperm(m-p)));
%
% where m is the length of the time series data U. Note that p observations
% are lost.
%
% REFS:
%
%    D. Freedman and D. Lane, "A Nonstochastic Interpretation of Reported
%    Significance Levels", J. Bus. & Econ. Stat. 1(4), pp. 292-298 (1983).
%
%    M.J. Anderson and J. Robinson, "Permutation tests for linear models",
%    Aust. N.Z. J. Stat. 43(1), pp. 75-88 (2001).

[n,m,N] = size(U);
assert (N == 1,'multi-trial version not yet implemented');

mp = m-p;
assert(mp > 0,'too many lags');
p1 = p+1;

U = demean(U);

z = 1:n; z([x y]) = [];    % omit y
conditional = ~isempty(z);
nx = length(x);
nxz = nx+length(z);

UU = U([x z y],:);         % rearrange

% full regression: X against lags of everything

X0 = UU(1:nx,p1:m);        % current observations
UL = zeros(p,n,mp);        % lagged  observations
for k = 1:p
    UL(k,:,:) = UU(:,p1-k:m-k); % k-lagged observations
end
UL  = reshape(UL,n*p,mp);  % stack lags
AX   = X0/UL;              % regress
AX = AX(:,1:p*nx);         % Axx
XX  = X0-AX*UL(1:p*nx,:);  % xi
AX = reshape(AX,nx,nx,p);  % in standard form

if conditional

    % reduced regression: xi against lags of Z

    ZL  = UL(p*nx+1:p*nxz,:);  % Z lags
    XZ = (XX/ZL)*ZL;           % predicted value of xi against Z
    EE  = XX-XZ;               % residuals eps'

else % unconditional

    % no need for reduced regression

    XZ = 0;
    EE = XX;

end
