function out = deterministic(in)
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
q = 1;

out.cidx = zeros(c,q);
out.ridx = zeros(r,q);

out.sigma_k = zeros(1,q);
out.froerr = zeros(1,q);
out.froerr_k = zeros(1,q);
out.specerr = zeros(1,q);
out.specerr_k = zeros(1,q);

out.construct_time = zeros(1,q);
out.metric_computing_time = zeros(1,q);

tic;
[~,~,Va]=svds(in.A,p);
out.cidx(:,1) = MSelect(Va(:,1:p),p,c);
C = in.A(:,out.cidx);

[~,~,Va]=svds(in.A',p);
out.ridx(:,1) = MSelect(Va(:,1:p),p,r);
R = in.A(out.ridx,:);
out.construct_time(1,1) = toc;

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

out.sigma_k(1,1) = Sb(end,end);
out.froerr(1,1) = norm(residual,'fro');
out.froerr_k(1,1) = norm(residual_k,'fro');
out.specerr(1,1) = svds(residual,1);
out.specerr_k(1,1) = svds(residual_k,1);

%out.trerr(1,iter) = trace(sqrt(residual*residual'));
%out.trerr(2,iter) = trace(sqrt(residual_k*residual_k'));

out.metric_computing_time(1,1) = toc;


end


