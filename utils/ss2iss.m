function [K,V,rep,L,P] = ss2iss(A,C,Q,R,S)

% Compute innovations form parameters for a state space model in general form by
% solution of a discrete algebraic Riccati equation (DARE) (eqs. 7, 8a, 8b in
% the reference article).
%
% A,C,Q,R,S - general form state space parameters
%
% K         - Kalman gain matrix
% V         - innovations covariance matrix
% rep       - DARE report (see below)
% L         - DARE stablising eigenvalues
% P         - DARE solution
%
% The value returned in rep is negative if an unrecoverable error was detected:
% rep = -1 means that the DARE was found to have eigenvalues on (or near) the
% unit circle, while rep = -2 indicates that no stabilising solution to the DARE
% could be found. If no error occurred, rep returns the relative residual of the
% DARE solution, which should be tested for accuracy (rep > sqrt(eps) is
% reasonable).
%
% Note that if the SS (A,C,Q,R,S) is well-formed - that is, A is stable and R
% positive definite - then (A,K,V) should be a well-formed innovations-form SS
% model. WE DON'T TEST FOR THAT HERE! It is up to the caller to do so if deemed
% necessary.

[r, r1] = size(A); assert(r1 == r);
[n, r1] = size(C); assert(r1 == r);
[r1,r2] = size(Q); assert(r1 == r && r2 == r);
if nargin < 5 || isempty(S)
	S = zeros(r,n);
else
	[r1,n1] = size(S); assert(r1 == r && n1 == n);
end
[n1,n2] = size(R); assert(n1 == n && n2 == n);
rr = 2*r;

K = [];
V = [];
P = [];

% We solve the DARE using the Generalized Schur (QZ) decomposition method

% Extended pencil

H = [A' zeros(r) C'; -Q  eye(r) -S; S' zeros(n,r) R];
J = [eye(r) zeros(r,r+n); zeros(r) A zeros(r,n); zeros(n,r) -C zeros(n)];

% NOTE - we don't balance!

[q,~] = qr(H(:,rr+1:rr+n));
H = q(:,n+1:rr+n)'*H(:,1:rr);
J = q(:,n+1:rr+n)'*J(:,1:rr);

% QZ algorithm

realHJ = isreal(H) && isreal(J);
i = 1:rr;
if realHJ
    [JJ,HH,q,z] = qz(J(i,i),H(i,i),'real');
else
    [JJ,HH,q,z] = qz(J(i,i),H(i,i),'complex');
end
[JJ,HH,~,z(i,:),qzflag] = ordqz(JJ,HH,q,z,'udo');
L = ordeig(JJ,HH);

% Check for stable invariant subspace

sis = abs(L) > 1;
if ~qzflag || any(~sis(1:r,:)) || any(sis(r+1:rr,:))
    rep = -1;
    return % IMPORTANT: caller must test!!! (error is: ''DARE: eigenvalues on/near unit circle'')
end

P1 = z(1:r,1:r);
P2 = z(r+1:rr,1:r);

% Solve P = P2/P1

[LL,UU,pvec] = lu(P1,'vector');
if rcond(UU) < eps
    rep = -2;
    return % IMPORTANT: caller must test!!! (error is: 'DARE: couldn''t find stabilising solution')
end
P(:,pvec) = (P2/UU)/LL;
P = (P+P')/2;

% Compute Kalman gain matrix K and innovations covariance matrix V

U = A*P*C'+S;
V = C*P*C'+R;
K = U/V;

% Check accuracy

APA = A*P*A'-P;
UK = U*K';
rep = norm(APA-UK+Q,1)/(1+norm(APA,1)+norm(UK,1)+norm(Q,1)); % relative residual

% IMPORTANT: test for accuracy  - something like
%
% if rep > sqrt(eps)
%     warning('DARE: possible inaccuracy (relative residual = %e)',rep);
% end

if nargout > 3

 % Return stable eigenvalues

    if realHJ
        L = L(r+1:rr);
    else
        L = conj(L(r+1:rr));
    end
end
