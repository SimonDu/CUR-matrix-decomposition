function [ idx ] = adaptive_sampling( A,idx21,r2 )
n = size(A,2);
R1 = A(idx21,:);
res = (A/R1)* R1;
res = A - res;
resNorm = zeros(n, 1);
for i = 1: n
    resNorm(i) = norm(res(:, i))^2;
end
clear res;
prob = resNorm / sum(resNorm);
idx22 = AdaptiveSampling(prob, r2);
idxtmp = 1: n;
idx21 = idxtmp(idx21);
idx = [idx21, idx22];
end

