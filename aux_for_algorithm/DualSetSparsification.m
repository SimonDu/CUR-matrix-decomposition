function [s] = DualSetSparsification(uu, V, r)
%%Dual Set Spectral-Frobenius Sparsification Algorithm
% input
%   U (m * n), V (k * n): two sets
%   r: number of non-zero entries of s
% output
%   s (n * 1): a set of weights containing at most r non-zero entries

[k n] = size(V);

% initialization
s = zeros(n, 1);
A = zeros(k, k);
sqrt_rk = sqrt(r*k);


j_old = 0;
for tau = 0: r-1
    Ltau = tau - sqrt_rk;
    [W, Lambda] = eig(A);    % W * Lambda * W' = A;
    lambda = diag(Lambda) - (Ltau + 1);
    dg_0 = 1 ./ (lambda + 1);
	dg_1 = 1 ./ lambda;
    clear Lambda lambda;
    dg_2 = dg_1 .^ 2;
    phi = sum(dg_1) - sum(dg_0);
        
    for idx = 0: n - 1
        j = mod(j_old + idx, n) + 1;
        v = W' * V(:, j);
        L = (v' * (dg_2 .* v)) / phi - ( v' * (dg_1 .* v) );
        
        clear v;
        
        if (uu(j) < L)
            t = 2 / (L + uu(j));
            s(j) = s(j) + t;
            A = A + t * (V(:, j) * V(:, j)');
            j_old = j;
            break;
        end
    end
    clear W dg_0 dg_1 dg_2;
end

%tmp = (1 - sqrt(k/r)) / r;
%s = tmp * s;


end