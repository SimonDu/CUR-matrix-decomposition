load abalone_dataset;
X = generate_distance_matrix(abaloneInputs');
in.A = generate_RBF_kernel(X, 5);
s5=svds(in.A,50);
in.A = generate_RBF_kernel(X, 2);
s2=svds(in.A,50);
in.A = generate_RBF_kernel(X, 0.2);
s02=svds(in.A,50);
in.A = generate_RBF_kernel(X, 0.1);
s01=svds(in.A,50);


