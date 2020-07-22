function plot_var(A,SIG,dt,scale,xlab,tstr)

[n,n1,p] = size(A);
assert(n1 == n,'VAR coefficients matrix has bad shape');
t = (1:p)'*dt;

if scale
    D = diag(sqrt(diag(SIG)));
    for k = 1:p
        A(:,:,k) = D\A(:,:,k)*D;
    end
end

xlims = [t(1) t(end)];

k = 0;
for i = 1:n
    for j = 1:n
        k = k+1;
        subplot(n,n,k);
        plot(t,squeeze(A(i,j,:)));
        grid on;
        if nargin > 2 && ~isempty(xlab), xlabel(xlab); end
        ylabel(sprintf('var %d,%d',i,j));
        xlim(xlims);
    end
end

if nargin > 5 && ~isempty(tstr), mtit(tstr); end
