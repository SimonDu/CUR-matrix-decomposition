function out = near_optimal(in)
%
% in is a structure with (at least) the following fields:
% - A, a matrix
% - k, the target rank of the approximation
% - c, number of columns to select
% - r, number of rows to select
% - q, the number of times to repeat each Nystrom method for each number of
%  column samples
% - sigma_k, 1 if we want output contain sigma_k, 0 otherwise
% - froerr, 1 if we want output contain froerr, 0 otherwise
% - froerr_k, 1 if we want output contain froerr_k, 0 otherwise
% - specerr, 1 if we want output contain specerr, 0 otherwise
% - specerr_k, 1 if we want output contain specerr_k, 0 otherwise
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
q = 1;
[m,n] = size(in.A);

out.cidx = {};
out.ridx = {};

if(in.sigma_k)
    out.sigma_k = zeros(1,q);
end
if(in.froerr)
    out.froerr = zeros(1,q);
end
if(in.froerr_k)
    out.froerr_k = zeros(1,q);
end
if(in.specerr)
    out.specerr = zeros(1,q);
end
if(in.specerr_k)
    out.specerr_k = zeros(1,q);
end

out.construct_time = zeros(1,q);
out.metric_computing_time = zeros(1,q);

for iter=1:in.q
    tic
    out.cidx{iter} = NearOptColSelect(in.A, in.k, c);
    C = in.A(:,out.cidx{iter});
    r1 = c;
    idx21 = NearOptColSelect(in.A', in.k, r1);
    R1 = in.A(idx21,:);
    res = (in.A/R1)* R1;
    r2 = r - r1;
    clear R1;
    if(r2 >0)
        res = in.A - res;
        resNorm = zero(n, 1);
        for i = 1: n
            resNorm(i) = norm(res(:, i))^2;
        end
        clear res;
        prob = resNorm / sum(resNorm);
        idx22 = AdaptiveSampling(prob, r2);
        idxtmp = 1: n;
        idx21 = idxtmp(idx21);
        out.ridx{iter} = [idx21, idx22];
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


