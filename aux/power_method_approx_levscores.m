function [approxlevscores, iter] = power_method_approx_levscores(A,k)
% [approxlevscores, iter] = power_method_approx_levscores(in)
%
% approximates the leverage scores of A, filtered through rank k,
% using the power method. The number of iterations is determined by testing
% that the approximations don't change by tol, or that a maximum number of
% iterations has been reached.
%
% in is a structure with (at least) the following fields:
% - A, an SPSD matrix
%
% chunk is the number of iterations after which to reorthogonalize and to check
% the leverage scores for convergence
%
% tol determines the tolerance for convergence: when the inf norm distance
% of the approximate leverage scores is smaller than this, convergence
% has been achieved
%
% maxiters is the maximum number of iterations to use
%
% returns the approximate leverage scores and the number of iterations
% needed
%
% Assumes k <= rank(A) and A is SPSD (so A' == A)

%these numbers are modified from gitten's code
maxiters = 10;
chunk = 10;
tol = 0.01;

n = size(A,1);
l = k;
S = randn(n,l);

Y = A*S;

iter = 0; % keep track of how many iterations of the power method were used

approxlevscores = zeros(1,n);

while iter < maxiters
    iter = iter + 1;
    if (rem(iter,chunk)==0)
        [Q,~] = qr(Y, 0);
        Y = Q; % With probability 1, Q still has k columns
        oldlevscores = approxlevscores;
        approxlevscores = sum(Q.^2,2)';
        if ( norm(approxlevscores - oldlevscores, Inf) < tol)
            return;
        end
    end
    Y = A*(A*Y); % since A is PSD, Y = A*A^T*Y
end

if (iter == maxiters)
    fprintf('Terminated with current approximate leverage scores after %i iterations\n', maxiters);
end
