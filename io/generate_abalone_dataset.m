function A = generate_abalone_dataset(sigma)
% A = generate_abalone_dataset(sigma)
%
% Takes: -sigma, RBF parameter
% Returns: -A, the similarity matrix, with A_{ij} =
% exp(-norm(x-y)^2/(2\sigma^2))
%
% Centers and normalizes the points before forming A
% Note that the average square distance between (normalized) points in the
% dataset is about 16

load abalone_dataset
numpts = size(abaloneInputs,2);
m = mean(abaloneInputs')';
stdv = std(abaloneInputs')';
pts = (abaloneInputs - repmat(m,1,numpts))./repmat(stdv,1,numpts);

A = zeros(numpts, numpts);
dists = A;
for row=1:size(pts,2)
    dists(row, :) = sum((pts - repmat(pts(:, row),1,numpts)).^2);
    A(row, :) = exp(-dists(row,:)/(2*sigma^2));
end

%disp(mean(dists(:)));
end