function out = uniform_sampling(in)
% out = simple_Nystrom(in)
%
% in is a structure with (at least) the following fields:
%  -A, an n-by-n PSD matrix, and
%  -k, the desired rank of the approximation
%  -l < n, the number of columns to sample and
%  -q >= 1 the number of times to repeat the experiment.
%
% out is a structure with the following fields:
%  - specerr, froerr, trerr: each matrices with two rows representing 
%  the spectral, frobenius, 
%  and trace norms of the errors in using q realizations of the simple 
%  column-based (sampled uniformly at random w/o replacement)
%  Nystrom extension to A, using l columns. The first row corresponds to 
%  Nystrom extensions where the rank was not fixed, the second to Nystrom
%  extensions where the rank was fixed
% 
%  -timings, a two row matrix of the time it took to run each experiment 
% (i.e. form C, Winv, Wkinv), including the time to approximate the leverage scores

% % if testing
%out = dummy_extension(in);
%return

n = size(in.A,1);
out.specerr = zeros(2,in.q);
out.froerr = zeros(2,in.q);
out.trerr = zeros(2,in.q);
out.timings = zeros(2,in.q);


for iter=1:in.q
    tic
        p = randperm(n);
        p = p(1:in.l);
        C = in.A(:, p);

        W = full(C(p,:));
    out.timings(:, iter) = toc;
    
    % time the eigenvalue decomposition of W
    tic
        [V,D] = orderedeig(W);
    decomptime = toc;
    
    % time to form Winv
    tic
        Winv = V*pinv(D)*V';
    out.timings(1, iter) = out.timings(1, iter) + decomptime + toc;
    
    % time to form Wkinv
    tic
        Wkinv = V(:, 1:in.k)*pinv(D(1:in.k, 1:in.k))*V(:,1:in.k)';
    out.timings(2, iter) = out.timings(2, iter) + decomptime + toc;
    
    [out.specerr(:,iter), out.froerr(:,iter), out.trerr(:,iter)] = ...
        estnorms(in, C, Winv, Wkinv);
end

end


