%setup
A = in.A;
[U,S,V] = svds(A);
target_singular_value = S(5,5);
fro_A_A_k = norm(A-U*S*V','fro');
spec_A_A_k = svds(A-U*S*V',1);
%p,c for plot
c = (50:10:100);