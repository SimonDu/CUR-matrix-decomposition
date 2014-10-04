function out = subspace_approxlevscore_gaussian(in)
%
% in is a structure with (at least) the following fields:
% - A, a matrix
% - k, the target rank of the approximation
% - p, the rank of first two partition matrix
% - c, number of columns to select
% - r, number of rows to select
% - q, the number of times to repeat each Nystrom method for each number of
%  column samples
%
% out is a structure with the following fields:
%  - cidx, c*q matrix represents the column index we choose for each
%  iteration
%  - ridx, r*q matrix represents the row index we choose for each
%  iteration
%  - sigma_k: 1*q vector represents the kth singular value of
%  reconstruction matrix for each iteration
%  - froerr: 1*q vector represents the error in frobenius norm of
%  reconstruction matrix for each iteration
%  - froerr_k: 1*q vector represents the error in frobenius norm of
%  truncated rank-k reconstruction matrix for each iteration
%  - specerr: 1*q vector represents the error in spectral norm of
%  reconstruction matrix for each iteration
%  - specerr_k: 1*q vector represents the error in spectral norm of
%  truncated rank-k reconstruction matrix for each iteration
%  - construct_time: 1*q vector represents the time to choose columns and rows
%  - metric_computing_time: 1*q vector represents the time to compute
%  different metrics

c = in.c;
r = in.r;
p = in.p;
q = in.q;

out.cidx = {};
out.ridx = {};

out.sigma_k = zeros(1,q);
out.froerr = zeros(1,q);
out.froerr_k = zeros(1,q);
% out.specerr = zeros(1,q);
% out.specerr_k = zeros(1,q);

out.construct_time = zeros(1,q);
out.metric_computing_time = zeros(1,q);
[m,n]=size(in.A);

for iter=1:in.q
    tic
    
    % compute the approximate leverage scores
    out.approxlevscores = ...
        spectral_sketch_levscores(in.A,p);
    levscoreprobs = out.approxlevscores/p;
    
    % sample according to those leverage scores
    colindices = ones(1,c);
    for i=1:c
        colindices(i) = find(cumsum(levscoreprobs) >= rand(),1);
    end
    scalingfactors = levscoreprobs(colindices).^(1/2);
    Id = speye(m);
    S = Id(:,colindices)*diag(1./scalingfactors);
    C = in.A*S;
    
    
    out.approxlevscores = ...
        spectral_sketch_levscores(in.A',p);
    levscoreprobs = out.approxlevscores/p;
    
    % sample according to those leverage scores
    colindices = ones(1,r);
    for i=1:r
        colindices(i) = find(cumsum(levscoreprobs) >= rand(),1);
    end
    scalingfactors = levscoreprobs(colindices).^(1/2);
    Id = speye(n);
    S = Id(:,colindices)*diag(1./scalingfactors);
    R = S'*in.A;
    
    out.construct_time = toc;
    
    tic
    [Qc,~] = qr(C,0);
    [Qr,~] = qr(R',0);
    
    B = Qc'*in.A*Qr;
    CUR = Qc*B*Qr';
    [Ub,Sb,Vb] = svds(B,in.k);
    Bk = Ub*Sb*Vb';
    CUR_k = Qc*Bk*Qr';
    
    residual = in.A-CUR;
    residual_k = in.A - CUR_k;
    
    out.specerr(1,iter) = svds(residual,1);
    out.specerr_k(1,iter) = svds(residual_k,1);
    out.froerr(1,iter) = norm(residual,'fro');
    out.froerr_k(1,iter) = norm(residual_k,'fro');
    out.sigma_k = Sb(end,end);
    out.metric_computing_time = toc;
    
end

end


