in.sigma = 2;
load abalone_dataset;
X = generate_distance_matrix(abaloneInputs');
in.A = generate_RBF_kernel(X, in.sigma);
clear X;

in.k = 20;
in.p = 39;
in.q = 15;

in.sigma_k = 1;
in.froerr = 1;
in.froerr_k = 1;
in.specerr = 0;
in.specerr_k = 0;

c_values = (in.p+1:10:150);
r_values = (in.p+1:10:150);


methods = {'deterministic','subspace_expected','subspace_approxlevscore_gaussian','subspace_approxlevscore_powermethod','uniform_sampling'};
out = run_different_number_of_c_and_r(in,methods,c_values,r_values);


save('./output/compare_abalone20_sigma_2')

c_values_plot;
export_fig(gcf,'./plots/c_plots_abalone_20_sigma_2.pdf');
close all;