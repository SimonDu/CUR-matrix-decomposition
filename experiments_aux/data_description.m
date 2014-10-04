function out = data_description(in)
% out = data_description(in)
% 
%
% A, a data matrix
% k, target rank
% p, oversampling parameter
%

[m,n] = size(in.A);
[U,S,V] = svds(in.A,max(in.p));
out.nonzero_percentage = nnz(in.A)/(m*n);


A_k = U(:,1:in.k)*S(1:in.k,1:in.k)*V(:,1:in.k)';
out.f_norm_ratio_k = norm(in.A-A_k,'fro')/norm(in.A,'fro');


out.sigma_p_sigma_k = [];
out.leverage_score_Vp = [];
out.coherence_score_Vp = [];
out.leverage_score_Up = [];
out.coherence_score_Up = [];


for i = 1:length(in.p)

   
    out.sigma_p_sigma_k(i) = S(in.p(i),in.p(i))/S(in.k,in.k);
    
    
    %out.n_norm_ratio_k = sum(svds(in.A-A_k,min(m,n)))/sum(diag(S));


    %calculate the leverage score of Vp
    L = 1:n;
    for j = 1:n
        L(j) = n/in.p(i)*(norm(V(j,1:in.p(i))))^2;
    end
    [sortedValues,~] = sort(L,'descend');
    out.leverage_score_Vp(i) = sortedValues(in.p(i));
    out.coherence_score_Vp(i) = max(L);
    %calculate the leverage score of Up
    L = 1:m;
    for j = 1:m
        L(j) = m/in.p(i)*(norm(U(j,1:in.p(i))))^2;
    end
    [sortedValues,~] = sort(L,'descend');
    out.leverage_score_Up(i) = sortedValues(in.p(i));
    out.coherence_score_Up(i) = max(L);
end
