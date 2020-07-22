% [H Z, ps, ps0, converged] = wilson_sf(S, fs, tol)
% Performs a numerical inner-outer factorization of a spectral matrix, using
% Wilsons method. This implementation here is a slight modification of the
% original implemention by M. Dhamala (mdhamala@gsu.edu) & G. Rangarajan
% (rangaraj@math.iisc.ernet.in) in August, 2006.
%
% modified by S K Mody (modysk@gmail.com), 22.Sept.2016
%
% Main reference used for the spectral Matrix Factorization: 
% The Factorization of Matricial Spectral Densities, SIAM J. Appl. Math,
% Vol. 23, No. 4, pgs 420-426 December 1972 by G T Wilson).
%
%% If you use this code for studying directed interactions (Granger causality), please cite
% the following references:
% -M.Dhamala, R.Rangarajan, M.Ding, Physical Review Letters 100, 018701 (2008)
% -M.Dhamala, R.rangarajan, M.Ding, Neuroimage 41, 354 (2008)
%
% ARGS:
% S:
%	Spectral matrix function. This should be specified as a (k x k x m)
%	array for the frequencies in the half open range [0, 0.5), divided
%	into equal intervals. The spectrum for negative frequencies is assumed
%	to be symmetric, ie:-
%		S(-f) = transpose(S(f))
%
% fs:
%	Sampling rate. This is required for normalization.
%	***IMPORTANT: Please ensure that the spectral matrix input, S, has
%	been normalized by the same value of fs, otherwise the the output
%	Z will be off by a factor.
%
% tol [default: 1e-6]:
%	The tolerance with which to check for convergence. Iterations stop
%	either when the number of iterations reaches a prespecified maximum
%	or when all of the following conditions are satisfied:-
%		|(ps0 - ps0_prev)./ps0_prev| < tol
%		|(ps - ps_prev)./ps_prev| < tol
%		|(S - ps*ps')./S| < tol
%	where |.| is the vector 2-norm.
%
% OUTPUT:
% H, Z:
%	H is complex array of the same size as S, and Z is real symmetric
%	positive definite matrix such that for each i:
%		S(:,:,i) = H(:,:,i)*Z*H(:,:,i)'
%
% ps, ps0:
%	(k x k x m) complex array. Theoretically Ps is a function defined on
%	the on the boundary of the unit circle in the complex plane, such that:
%		S(:,:,i) = ps(:,:,i)*ps(:,:,i)'
%	Theoretically, Ps has a holomorphic extension in the complex plane to
%	all |z| < 1.). ps0 is the upper triangular matrix that is the value of
%	ps at the origin. Z is related to ps0 by:
%		Z = ps0*ps0'
%	
% converged:
%	Boolean value indicating whether the iteration converged to within the
%	specified tolerance.
%
% relerr:
%	The relative Cauchy error of the convergence of the spectrum or Ps.
%
function [H, Z, ps, ps0, converged, relerr] = wilson_sf(S, fs, tol)
	if (nargin < 3) || isempty(tol), tol = 1e-6; end
	assert(isscalar(fs) && (fs > 0), ...
		'fs must be a positive scalar value representing the sampling rate. ');
	
	[k, ~, N] = size(S);

	Sarr = cat(3, S, conj(S(:, :, N:-1:2)));
	ps0 = ps0_initial__(Sarr);
	ps = repmat(ps0, [1,1,N]);
	ps = cat(3, ps, conj(ps(:,:,N:-1:2)));
	M = 2*N-1;

	I = eye(k);
	maxiter = min( 500, floor(sqrt(10/tol)) );

	U = zeros(size(Sarr));
	for j = 1 : M
		U(:,:,j) = chol(Sarr(:,:,j));
	end
	%V = zeros(size(Sarr));

	niter = 0;
	converged = false;
	g = zeros(k,k,M);
	while ( (niter < maxiter) && ~converged )
		for i = 1 : M
			% Equivalent to:
			% g(:,:,i) = ps(:,:,i)\Sarr(:,:,i)/ps(:,:,i)' + I;
			V = ps(:,:,i)\U(:,:,i)';
			g(:,:,i) = V*V' + I;
		end

		[gp, gp0] = PlusOperator(g);
		T = -tril(gp0, -1);
		T = T - T';

		ps_prev = ps;
		for i = 1 : M,
			ps(:,:,i) = ps(:,:,i)*(gp(:,:,i) + T);
		end

		ps0_prev = ps0;
		ps0 = ps0*(gp0 + T);

		% Relative cauchy error. Check on S is expensive, so check Ps0 first, then Ps and only then S.
		[converged relerr] = check_converged_ps__(ps, ps_prev, ps0, ps0_prev, tol);
		if converged
			% Uncomment this next line to check for relative cauchy error in spectrum.
			%[converged relerr] = check_converged_S__(Sarr, ps, tol);
		end

		niter = niter + 1;
	end

	H = zeros(k,k,N);
	for i = 1 : N,
		H(:,:,i) = ps(:,:,i)/ps0;
	end
	
	ps = sqrt(fs)*ps(:,:,1:N);
	ps0 = sqrt(fs)*ps0;
	Z = ps0*ps0';
end

function ps0 = ps0_initial__(Sarr)
	[k, ~, M] = size(Sarr);
	
	% perform ifft to obtain gammas.
	Sarr = reshape(Sarr, [k*k, M]);
	gamma = ifft(transpose(Sarr));
	gamma0 = gamma(1,:);
	gamma0 = reshape(gamma0, [k k]);
	
	% Remove any assymetry due to rounding error.
	% This also will zero out any imaginary values
	% on the diagonal - real diagonals are required for cholesky.
	gamma0 = real((gamma0 + gamma0')/2);

	ps0 = chol(gamma0);
	
end

%% This function is for [ ]+operation
function [gp, gp0] = PlusOperator(g)

	[k, ~, M] = size(g);
	N = (M+1)/2;
	
	g = reshape(g, [k*k, M]);
	gammma = real(ifft(transpose(g)));
	gammma = reshape(transpose(gammma), [k,k,M]);
	
	% Take half of the zero lag
	gammma(:,:,1) = 0.5*gammma(:,:,1);
	gp0 = gammma(:,:,1);
	
	% Zero out negative powers.
	gammma(:, :, N+1:end) = 0;

	% Reconstitute
	gammma = reshape(gammma, [k*k, M]);
	gp = fft(transpose(gammma));
	gp = reshape(transpose(gp), [k,k,M]);
	
end
%%

function [converged_ps relerr] = check_converged_ps__(ps, ps_prev, ps0, ps0_prev, tol)

	[converged_ps relerr] = CheckRelErr__(ps0, ps0_prev, tol);
	if converged_ps
		[converged_ps RelErr2] = CheckRelErr__(ps, ps_prev, tol);
		relerr = max(relerr, RelErr2);
	end
	
end
%%

function [converged_S relerr] = check_converged_S__(S, ps, tol)

	FX = zeros(size(ps));
	parfor j = 1 : size(ps,3)
		FX(:,:,j) = ps(:,:,j)*ps(:,:,j)';
	end
	[converged_S relerr] = CheckRelErr__(FX, S, tol);
	
end
%%

function [ok, relerr] = CheckRelErr__(A,B,reltol)
	d = norm(B(:)-A(:), 2);
	a = norm(A(:), 2);
	
	if a > 0, relerr = d/a;
	else relerr = d;
	end
	
	ok = (relerr <= reltol);
end
%%

