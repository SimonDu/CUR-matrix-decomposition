function [V] = GaussProjSVDRight(A, k, p)
%%Compute truncated SVD using Gaussian projection
% Input:
%   A: data matrix at hand
%   k: target rank
%   p: oversampling parameter (default p = k)
% Output:
%   V: top k right singular vectors of A

if nargin < 3
    l = 4 * k;
else
    l = k + p;
end

[n] = size(A, 2);

%% Sate A
% Gaussian projection matrix
Omega = randn(n, l);

% Columns of Q form the range of Y
Y = A * Omega;
[Q,~] = qr(Y,0);
%[Q S V] = svd(Y, 'econ');

%% Stage B
B = Q' * A;
[~, ~,  V] = svd(B, 'econ');

V = V(:, 1:k);
end