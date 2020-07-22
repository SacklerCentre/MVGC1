function plot_infocrit(mos,AIC,BIC,dt,xlab,tstr)

if isscalar(mos), mos = 1:mos; end
assert(isvector(mos),'model orders must be a scalar or a vector');
nmos = length(mos);
assert(isvector(AIC) && isvector(BIC) && length(AIC) == nmos && length(BIC) == nmos,'AIC and BIC must be vectors matching the model orders');

t = mos(:)*dt;
xlims = [t(1) t(end)];

plot(t,[AIC(:) BIC(:)]);
xlim(xlims);
legend('AIC','BIC');
if nargin > 4 && ~isempty(xlab), xlabel(xlab); end
if nargin > 5 && ~isempty(tstr), title(tstr); end
