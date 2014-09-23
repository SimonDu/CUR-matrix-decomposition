function out = near_optimal(in)
%
% in is a structure with (at least) the following fields:
% - A, a matrix
% - k, the target rank of the approximation
% - c, number of columns to select
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
q = in.q;
[~,n] = size(in.A);

out.cidx = {};
out.ridx = {};

out.sigma_k = zeros(1,q);
out.froerr = zeros(1,q);
out.froerr_k = zeros(1,q);
out.specerr = zeros(1,q);
out.specerr_k = zeros(1,q);

out.construct_time = zeros(1,q);
out.metric_computing_time = zeros(1,q);
for iter=1:in.q
    tic
        idx1 = NearOptColSelect(in.A, in.k, c);
        C = in.A(:,idx1);
        epsilon = 2*in.k/c;
        r = c;
        idx21 = NearOptColSelect(in.A', in.k, r);
        R1 = in.A(idx21,:);
        res = (in.A/R1)* R1;
        r2 = c/epsilon;
        clear R1;
        res = in.A - res;

        resNorm = ones(n, 1);
        for i = 1: n
            resNorm(i) = norm(res(:, i))^2;
        end
        clear res;
        prob = resNorm / sum(resNorm);

        idx22 = AdaptiveSampling(prob, r2);

        idxtmp = 1: n;
        idx21 = idxtmp(idx21);
        idx2 = [idx21, idx22];
        R = in.A(idx2,:);
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
        %out.trerr(1,iter) = trace(sqrt(residual*residual'));
        %out.trerr(2,iter) = trace(sqrt(residual_k*residual_k'));
        out.sigma_k = Sb(end,end);
    out.metric_computing_time = toc;
    
end

end


