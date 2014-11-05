function out = randomized_unweighted(in)
% Randomized unweighted column selection based algorithm for CUR
% decomposition
% in is a structure with (at least) the following fields:
% - A, a matrix
% - k, the target rank of the approximation
% - p, the rank of first two partition matrix
% - c, number of columns to select
% - r, number of rows to select
% - adaptive, 1 if we want to do adaptive sampling(assume r > c), 0
% othewise(r = c)
% - sigma_k, 1 if we want output contain sigma_k, 0 otherwise
% - froerr, 1 if we want output contain froerr, 0 otherwise
% - froerr_k, 1 if we want output contain froerr_k, 0 otherwise
% - specerr, 1 if we want output contain specerr, 0 otherwise
% - specerr_k, 1 if we want output contain specerr_k, 0 otherwise

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

k = in.k;
c = in.c;
r = in.r;
p = in.p;

out.cidx = {};
out.ridx = {};

if(in.sigma_k)
    out.sigma_k = zeros(1,1);
end
if(in.froerr)
    out.froerr = zeros(1,1);
end
if(in.froerr_k)
    out.froerr_k = zeros(1,1);
end
if(in.specerr)
    out.specerr = zeros(1,1);
end
if(in.specerr_k)
    out.specerr_k = zeros(1,1);
end

out.construct_time = zeros(1,1);
out.metric_computing_time = zeros(1,1);
[m,n] = size(in.A);

for iter = 1:1
    tic;
    Omega = randn(n,2*p);
    Y = in.A*Omega;
    for i = 1:1
        Y = in.A*(in.A'*Y);
    end
    [Q,~] = qr(Y,0);
    B = Q'*in.A;
    [~,~,Va] = svd(B,'econ');
    out.cidx{iter} = MSelect(Va(:,1:p),p,c);
    C = in.A(:,out.cidx{iter});
    
    Omega = randn(m,2*p);
    AT = in.A';
    Y = AT*Omega;
    for i = 1:1
        Y = AT*(AT'*Y);
    end
    [Q,~] = qr(Y,0);
    B = Q'*AT;
    [~,~,Va]=svd(B,'econ');
    idx21 = MSelect(Va(:,1:p),p,c);
    if(in.adaptive)
        out.ridx{iter} = adaptive_sampling(in.A,idx21,r-c);
    else
        out.ridx{iter} = idx21;
    end
    R = in.A(out.ridx{iter},:);
    
    out.construct_time(1,iter) = toc;
    
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
    
    if(in.sigma_k)
        out.sigma_k(1,iter) = Sb(end,end);
    end
    if(in.froerr)
        out.froerr(1,iter) = norm(residual,'fro');
    end
    if(in.froerr_k)
        out.froerr_k(1,iter) = norm(residual_k,'fro');
    end
    if(in.specerr)
        out.specerr(1,iter) = svds(residual,1);
    end
    if(in.specerr_k)
        out.specerr_k(1,iter) = svds(residual_k,1);
    end
    
    out.metric_computing_time(1,iter) = toc;
end


end


