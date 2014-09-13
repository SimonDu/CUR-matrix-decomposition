function [sparseA, L, nL] = read_snap_data(fname)
% [sparseA, L, nL] = read_snap_data(fname)
% 
% Reads in (some) files from the Stanford Large Network Dataset collection
% (SNAP) and returns the corresponding adjacency matrix, unnormalized
% Laplacian, and normalized Laplacian
%
% Takes:
% -fname, the name of a file containing undirected (but symmetric) data
% from the Stanford Large Network Dataset collection
%
% Returns:
% -sparseA, the adjacency matrix of the graph, in sparse matrix form, and
% -L, the sparse matrix unnormalized Laplacian of the graph, and
% -nL, the sparse matrix normalized Laplacian of the graph (I - D^{-1/2}A
% D^{-1/2})
%

initskip = 4;
edgepairs = importdata(fname,'\t',initskip);

% the nodes aren't numbered 1 through numnodes, so fix this
% [b,~,j] = unique(edgepairs.data(:,1));
% [c,~,k] = unique(edgepairs.data(:,2));
% n = length(b); % number of nodes
% assert(norm(b-c)==0);
% noderange = [1:n]';
% fixededgepairs = [noderange(j), noderange(k)];

[orignodelabels,~,nodeordering] = ...
    unique([edgepairs.data(:,1); edgepairs.data(:,2)]);
n = length(orignodelabels);
numpairs = size(edgepairs.data,1);
noderange = [1:n]';
fixededgepairs = [noderange(nodeordering(1:numpairs)), ...
    noderange(nodeordering((numpairs+1):end))];

% The header in the SNAP files say that there's only one directed edge entry per
% unordered edge, so we should have to take sparseA = sparseA + sparseA' to make it an
% adjacency matrix, but in fact they record both directed edges per
% unordered edge, so it's fine to just take sparseA as follows
sparseA = sparse(fixededgepairs(:,1), fixededgepairs(:,2),1);

% Added later (for p2p-gnutella dataset, sparseA is actually not symmetric)
if svds(sparseA - sparseA',1) > 10*eps
    fprintf('Forcing symmetry in adjacency graph\n');
    sparseA = sparseA + sparseA';
end

% form unnormalized and normalized Laplacians
degrees = sum(sparseA);
Dnegsqrt = diag(degrees.^(-1/2));
L = sparse(diag(degrees)) - sparseA;
nL = sparse(eye(n)) - Dnegsqrt*sparseA*Dnegsqrt;

end
