[m,n] = size(in.A)
nnz_ratio = 100*nnz(in.A)/(m*n)
k = in.k
[U,S,V] = svds(in.A,k);
A_k = U*S*V';


f_norm = norm(in.A,'fro');
two_norm = svds(in.A,1);

f_2_ratio = (f_norm/two_norm)^2
A_k_left_ratio = 100*norm(in.A-A_k,'fro')/norm(in.A,'fro')

for i=1:m
    L(i) = norm(U(i,:));
end
mu_c = m/k * max(L)

for i=1:n
    L(i) = norm(V(i,:));
end
mu_r = n/k * max(L)

s = svds(in.A,2*k);
sigma_k_ratio = s(2*k)/s(k)