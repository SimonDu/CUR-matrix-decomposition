function out = deterministic(in)
% out = deterministic(in)
%
% in is a structure with (at least) the following fields:
% - A, a matrix
% - k, the target rank of the approximation
% - p, the rank of first two partition matrix
% - number_of_c_and_r, a two row matrix specifying the numbers of column and rows 
%   samples to use   
% - q, the number of times to repeat each Nystrom method for each number of
%  column samples
%
% out is a structure with the following fields:
%  - cidx, index vector of selected columns
%  - ridx, index vector of selected rows
%  - specerr, froerr, trerr: each matrices with two rows representing 
%  the spectral, frobenius, 
%  and trace norms of the errors in using q realizations of the simple 
%  column-based (sampled uniformly at random w/o replacement)
%  Nystrom extension to A, using l columns. The first row corresponds to 
%  Nystrom extensions where the rank was not fixed, the second to Nystrom
%  extensions where the rank was fixed
%  - sigma_k, the k-th singular value of A
% 
%  -timings, a two row matrix of the time it took to run each experiment 
% (i.e. form C, Winv, Wkinv), including the time to approximate the leverage scores

% % if testing
%out = dummy_extension(in);
%return


[m,n] = size(in.A);
c = in.number_of_c_and_r(1);
r = in.number_of_c_and_r(2);
out.specerr = zeros(2,in.q);
out.froerr = zeros(2,in.q);
out.trerr = zeros(2,in.q);
out.timings = zeros(2,in.q);

for iter=1:in.q
    tic
        [~,~,Va]=svds(in.A,in.p);
        cidx = MSelect(Va(:,1:in.k),in.k,c);
        C = in.A(:,cidx);

        [~,~,Va]=svds(in.A',in.p);
        ridx = MSelect(Va(:,1:in.k),in.k,r);
        R = in.A(ridx,:);
    out.timings(1, iter) = toc;
    
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
        out.specerr(2,iter) = svds(residual_k,1);
        out.froerr(1,iter) = norm(residual,'fro');
        out.froerr(2,iter) = norm(residual_k,'fro');
        %out.trerr(1,iter) = trace(sqrt(residual*residual'));
        %out.trerr(2,iter) = trace(sqrt(residual_k*residual_k'));
        out.sigma_k = Sb(end,end);
    out.timings(2,iter) = toc;

end

end


