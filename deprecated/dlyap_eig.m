function G = dlyap_eig(A,SIG)

n = size(A,1);
[U,D] = eig(A);
V = U\eye(n);
d = diag(D);
G = real(U*((V*SIG*V')./(1-d*d'))*U');
