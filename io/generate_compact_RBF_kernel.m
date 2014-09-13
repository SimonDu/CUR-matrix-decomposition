function A = generate_compact_RBF_kernel(distmat, sigma, d, cutoff)
% A = generate_compact_RBF_kernel(distmat, sigma, d, cutoff)
%
% Returns: 
% -A, the similarity matrix, with A_{ij} = c(x,y,cutoff)*exp(-norm(x-y)^2/sigma^2)
%
% Takes: 
% -distmat, (ij)th entry is norm(x_i - x_j)^2
% -sigma, RBF parameter
% -d, the dimension of the points x_i
% -cutoff, a sparsification parameter
%
% where c(x,y,cutoff) = max(0, (1 - norm(x-y)/cutoff)^ceil((d+1)/2) )
%
% A standard choice of cutoff is e.g. 3*sigma
%
% Reference: Genton 2001. Classes of Kernels for Machine Learning: A Statistics
% Perspective

A = zeros(size(distmat));
p = ceil((d+1)/2);
c = @(rowofdists) max([1 - sqrt(rowofdists)/cutoff; zeros(size(rowofdists))]).^p;

for row=1:size(distmat,1)
    A(row, :) = c(distmat(row,:)).*exp(-distmat(row,:)/sigma^2);
end

end
