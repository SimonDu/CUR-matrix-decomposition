function Xhat = normalize_kernel_data(X)
% function Xhat = normalize_kernel_data(X)
%
% given an n-by-d matrix of n observations in R^d,
% normalizes X so that A = Xhat*Xhat' is a linear kernel matrix with
% nice properties. Returns Xhat.
%
% namely: it first centers the columns of X and sets them to have std dev 1
% so all the observations on the same scale of importance
% then it scales the rows of X so they all have norm 1 --- this
% ensures the diagonals of the kernel are 1

numpts = size(X,1);
numfeatures = size(X,2);

% make observations zero-mean and variance one
m = full(mean(X));
stdv = full(std(X));
stdv(stdv < 10^(-4)) = 1; % avoid division by 0
Xhat = (X - repmat(m, numpts, 1))./repmat(stdv, numpts, 1);

% normalize the diagonals of the kernel
Xhat = Xhat ./ repmat(sqrt(sum(Xhat.*Xhat, 2)), 1, numfeatures);

end
