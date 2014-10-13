in.sigma = 2;
load abalone_dataset;
X = generate_distance_matrix(abaloneInputs');
in.A = generate_RBF_kernel(X, in.sigma);
clear X;

in.k = 20;
in.p = 39;
in.q = 10;

in.sigma_k = 1;
in.froerr = 1;
in.froerr_k = 1;
in.specerr = 0;
in.specerr_k = 0;

in.adaptive = 1;

c_values = (in.p+1:20:140);
r_values = 2*c_values;


methods = {'deterministic','randomized_unweighted','subspace_expected','subspace_approxlevscore_gaussian','uniform_sampling','near_optimal'};
out = run_different_number_of_c_and_r(in,methods,c_values,r_values);


save('./output/compare_abalone20_sigma_2_adaptive')

c_values_plot;
saveas(gcf,'./plots/c_plots_abalone_20_sigma_2_adaptive','fig');
export_fig(gcf,'./plots/c_plots_abalone_20_sigma_2_adaptive.pdf');
close all;