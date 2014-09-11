function [snorm, fnorm, tnorm] = estnorms(in, C, Winv, Wkinv)
% [snorm, fnorm, tnorm] = estnorms(in, C, Winv, Wkinv)
%
% Given A, k, C, and Winv, Wkinv, estimates the spectral, frobenius, and trace norms of 
% A - C*Winv*C' and A - C*Wkinv*C'
% where Winv is the pseudoinverse of W and Wkinv is the pseudoinverse of an
% optimal rank k approximation to W.
%
% in is a structure with (at least) the following fields:
% A -- an SPSD matrix
% k -- the rank of the fixed-rank approximation
%
% snorm, fnorm, tnorm are vectors whose first entries represent the errors of
% using Winv and whose second entries represent the errors of using Wkinv

% NOTE: normest is horrible (it says the norm of the enron laplacian matrix
%  is 1.98... when it is actually 2), and I've come across cases where eigs 
% says a very large nonzero matrix has norm zero, so switch to using 
% jdrpcg (jdqr is more appropriate, but it has all sorts of implementation
% issues: not returning all the eigenvalues you ask for, failing on some
% weird cases...)
% so use a try catch to use eigs then fall back to jdrpcg if that fails

% first calculate the errors of the non-fixed rank approximation

D = in.A - C*Winv*C';

try
    snorm(1) = eigs(D, 1);
    if snorm(1) < 0
        error('The input matrix is not positive!');
    end
catch err
    try
    	fprintf('Warning: switching to using jdrpcg because eigs failed to converge\n');
    	snorm(1) = -jdrpcg(-D, 1);
    	if snorm(1) < 0
        	error('The input matrix is not positive!');
    	end
    catch err
	fprintf('Truly exceptional case!: jdrpcg failed, so using norm\n');
	snorm(1) = norm(D);
    end
end

fnorm(1) = norm(D, 'fro');
tnorm(1) = trace(D); %appropriate since D is presumed SPSD

% next calculate the errors of the fixed rank approximation

D = in.A - C*Wkinv*C';

%snorm(2) = jdqr(D,1); % appropriate since A is presumed SPSD
try
	snorm(2) = eigs(D,1);
    if snorm(2) < 0
        error('The input matrix is not positive!');
    end
catch err
	snorm(2) = -jdrpcg(-D,1);
    fprintf('Warning: switching to using jdrpcg because eigs failed to converge\n');
    if snorm(2) < 0
        error('The input matrix is not positive!');
    end
end
fnorm(2) = norm(D, 'fro');
tnorm(2) = trace(D); %appropriate since D is presumed SPSD

snorm = real(snorm); % sometimes have spurious 0i components
return;

