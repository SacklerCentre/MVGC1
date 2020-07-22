function [tstat,gamma] = adf(x,pdeg,q)

%---------------------------------------------------
% FUNCTION: cca_adf(x,pdeg,q)
% INPUTS:   x = a time-series vector
%           pdeg = order of time polynomial in the null-hypothesis
%           pdeg = -1, no deterministic part
%           pdeg =  0, for constant term
%           pdeg =  1, for constant plus time-trend
%           pdeg >  1, for higher order polynomial
%           nlags = # of lagged changes of x included
%
% OUTPUT: a results structure
%         results.meth  = 'adf'
%         results.gamma = estimate of the autoregressive parameter
%         results.adf   = ADF t-statistic
%         results.crit = (6 x 1) vector of critical values
%                        [1% 5% 10% 90% 95% 99%] quintiles
%---------------------------------------------------
% References: Said and Dickey (1984) 'Testing for Unit Roots in
% Autoregressive Moving Average Models of Unknown Order',
% Biometrika, Volume 71, pp. 599-607.
%
% written by:
% James P. LeSage
% Modeled after a similar Gauss routine by
% Sam Ouliaris, in a package called COINT
%
% modified and encapsulated by Anil Seth, December 2005
%---------------------------------------------------

assert(pdeg >= -1 && pdeg <= 5,'trend parameter must lie between -1 and 5');
assert(isvector(x),'input time series must be a vector');

x = x(:); % make sure it's a column vector;
m = length(x);
assert(2*q <= m,'too many lags');
q1 = q+1;

% Gerard van den Hout suggested the fix below
% Erasmus University Rotterdam.
% The Netherlands.

DX = diff(x,1);
z = zeros(m-1,q);
for k = 1:q
    z(:,k) = [zeros(k,1); DX(1:end-k,:)];
end
z = z(q1:end,:);
if pdeg > -1
    z = [z ptrend(pdeg,m-q1)];
end
zcov = z'*z;
nz = size(zcov,1);
izcov = zcov\eye(nz);
dep = x(q+2:end,:);
b   = izcov*z'*dep;

% res     = dep - z*b ;
% BUG fix suggested by
% Nick Firoozye
% Sanford C. Bernstein, Inc

res  = tdemean(dep) - tdemean(z)*b;
so   = (res'*res)/(m-q1-nz);

gamma = b(1,1);
tstat = (b(1,1)-1)/sqrt(so*izcov(1,1));

%----------------------------------------------------------
function xmat = ptrend(pdeg,m)

    if pdeg > 0
        timep = zeros(m,pdeg) ;
        t = (1:m)'/m;
        for r = 1:pdeg
            timep(:,r) = t.^r ;
        end
        xmat = [ones(m,1) timep];
    else
        xmat = ones(m,1);
    end

%----------------------------------------------------------
function z = tdemean(y) % temporal demean

    z = y-ones(size(y,1),1)*mean(y);
