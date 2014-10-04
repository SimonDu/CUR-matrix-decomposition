function approxlevscores = spectral_sketch_levscores(A,k)

% approxlevscores = spectral_sketch_levscores(in)
%
% approximates the leverage scores of A, filtered through rank k,
% using Algorithm 4 of Mahoney et al.
%
% in is a structure with (at least) the following fields:
% - A
% - k
%
%
% Assumes k <= rank(A)


% chunk is the number of iterations after which to reorthogonalize
% returns the approximate leverage scores and the number of iterations
% needed

chunk = 10;
n = size(A,1);
S = randn(n, 2*k);
% Estimate 2*log(1 + eps/10) - 1/2 from the bound in the paper as 1
% q = ceil( log(1 + sqrt(k/(k-1)) + exp(1)*sqrt(2*(n-k)/k)) );
q = 4; % if you use the above expression, for n=10000, k=30, q would be about 5



Y = A*S;

iter = 0; % keep track of how many iterations of the power method were used

while iter < q
    iter = iter + 1;
    % reorthogonalize every chunk steps
    if (rem(iter,chunk)==0)
        [Q,~] = qr(Y, 0);
        Y = Q; % With probability 1, Q still has l columns
    end
    Y = A*(A'*Y); % since A is PSD, Y = A*A^T*Y
end

[U,~,~] = svds(Y, k);
approxlevscores = sum(U.^2, 2)';