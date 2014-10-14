[m,n] = size(in.A)
nnz(in.A)/(m*n)
k = 10
[U,S,V] = svds(in.A,k);
for i=1:m
    L(i) = norm(U(i,:));
end
mu_c = m/k * max(L)

for i=1:n
    L(i) = norm(V(i,:));
end
mu_r = n/k * max(L)