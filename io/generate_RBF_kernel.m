function A = generate_RBF_kernel(distmat, sigma)
% A = generate_RBF_kernel(distmat sigma)
%
% Returns: 
% -A, the similarity matrix, with A_{ij} = exp(-norm(x-y)^2/sigma^2)
%
% Takes: 
% -distmat, (ij)th entry is norm(x_i - x_j)^2
% -sigma, RBF parameter

A = zeros(size(distmat));
for row=1:size(distmat,1)
    A(row, :) = exp(-distmat(row,:)/sigma^2);
end

end
