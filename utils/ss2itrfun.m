function J = ss2itrfun(A,C,K,fres)

[n,r] = size(C);
h = fres+1;
In = eye(n);
Ir = eye(r);
B = A-K*C;
J = zeros(n,n,h);
w = exp(1i*pi*((0:fres)/fres));
for k = 1:h % over [0,pi]
    J(:,:,k) = In - C*((w(k)*Ir-B)\K); % NOTE: if minimum phase then w(k)*Ir-B always invertible!
end
