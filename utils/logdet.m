%% logdet
%
% Calculate logarithm of determinant for positive-definite matrices
%
% <matlab:open('logdet.m') code>
%
%% Syntax
%
%     L = logdet(V)
%
%% Arguments
%
% _input_
%
%     V          a positive-definite square, symmetric matrix
%
% _output_
%
%     L          log-determinant of V
%     flag       zero if V is Hermitian positive-definite, else nonzero
%
%% Description
%
% Returns the log-determinant of the Hermitia positive-definite matrix V in L.
%
% Essentially the same as |log(det(V))|, but avoids potential under/overflow. The 'flag'
% variable is zero if V is Hermitian positive-definite. If V is not Hermitian positive-
% definite, the real part of the log-determinant is returned if the imaginary part is not
% too large, else NaN.
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function [LD,flag] = logdet(V)

[L,flag] = chol(V);
if flag == 0 % symmetric, positive-definite
	LD = 2*sum(log(diag(L)));
else
	DV = det(V);
	if abs(imag(DV)) > sqrt(eps)
		LD = NaN;
	else     % give it the benefit of the doubt...
		LD = log(real(DV));
	end
end
