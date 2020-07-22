function plot_morder(mos,pval,dt,amax,ares,tstr)

if isscalar(mos), mos = 1:mos; end
assert(isvector(mos),'model orders must be a scalar or a vector');
nmos = length(mos);
assert(isvector(pval) && length(pval) == nmos,'p-values must be a vector matching the model orders');

if nargin < 4 || isempty(amax), amax = 0.1; end
if nargin < 5 || isempty(ares), ares = 100; end

subplot(2,1,1);
t = mos(:)*dt;
plot(t,pval);
xlim([t(1) t(end)]);
ylim([0,2*amax]);
xlabel('model order');
ylabel('p-value');
grid on

alpha = linspace(0,amax,ares);
moa = nan(ares,1);
for j = 1:ares
    [~,moa(j)] = min(pval < alpha(j));
end
subplot(2,1,2);
plot(alpha,moa);
xlabel('significance level');
ylabel('model order');
grid on

if nargin > 5 && ~isempty(tstr), mtit(tstr); end
