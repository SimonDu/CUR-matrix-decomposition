[m,n] = size(in.A)
100*nnz(in.A)/(m*n)
k = 10
[U,S,V] = svds(in.A,k);
A_k = U*S*V';

f_norm = norm(in.A,'fro');
two_norm = svds(in.A,1);

ceil(f_norm/two_norm)
100*norm(in.A-A_k,'fro')/norm(in.A,'fro')

for i=1:m
    L(i) = norm(U(i,:));
end
mu_c = m/k * max(L)

for i=1:n
    L(i) = norm(V(i,:));
end
mu_r = n/k * max(L)

