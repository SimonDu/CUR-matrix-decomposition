in.sigma = 1;
X = generate_abalone_dataset(in.sigma); 
in.A = generate_RBF_kernel(X, in.sigma);
clear X;

in.k = 20;
in.c = floor(2*in.k*log(in.k));
in.r = floor(2*in.k*log(in.k));
in.q = 10;

methods = {'subspace_expected','deterministic'};
p_values = (20:80);

out = run_dataset_different_p(in,methods,p_values);

p_values_plot_deterministic;
p_values_plot_subspace_expected;