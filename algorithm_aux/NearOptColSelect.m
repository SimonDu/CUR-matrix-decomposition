function [idx] = NearOptColSelect(A, k, c)
%%Near Optimal Column Selection
% input
%   A (m*n): the given matrix to decompose
%   k:  target rank
%   c:  number of columns to be selected
% output
%   idx: selected columns

[m n] = size(A);


%% SVD via Random Projection
% over-sampling parameter
l = 4*k;
if l >= min(m, n)
    [UA SA VA] = svds(A, k);
    clear SA UA;
else
    [VA] = GaussProjSVDRight(A, k, l-k);
end

Ares = A * VA * VA';
Ares = A - Ares;



%% Set parameters
epsilon = 2 * k / c;
c1 = 2 * k * epsilon^(-2/3);
c1 = floor(c1);
c = min(c1, c-1);
clear epsilon;


%% Dual Set Sparsifications
AresFroSq = norm(Ares, 'fro')^2;
deltaU1 = AresFroSq / (1 - sqrt(k / c));

uu1 = zeros(n, 1);
for iU = 1: n
    uu1(iU) = norm(Ares(:, iU))^2;
end
uu1 = uu1 / deltaU1;
clear Ares deltaU1 AresFroSq iU;

s = DualSetSparsification(uu1, VA', c);
clear VA uu1;
idx11 = (s' > 0);


%% Adaptive Sampling
C1 = A(:, idx11);
res = C1 * (C1 \ A);

clear C1;
res = A - res;


resNorm = ones(n, 1);
for i = 1: n
    resNorm(i) = norm(res(:, i))^2;
end
clear res;
prob = resNorm / sum(resNorm);

c2 = c - sum(idx11);
idx12 = AdaptiveSampling(prob, c2);

idxtmp = 1: n;
idx11 = idxtmp(idx11);
idx = [idx11, idx12];
end


