% RBF kernel for abalone dataset with sigma = 5 and target rank = 20
in.sigma = 5;
load abalone_dataset;
X = generate_distance_matrix(abaloneInputs');
in.A = generate_RBF_kernel(X, in.sigma);
clear X;

in.k = 20;
in.p = 30;
in.q = 10;

in.sigma_k = 1;
in.froerr = 1;
in.froerr_k = 1;
in.specerr = 1;
in.specerr_k = 1;

in.adaptive = 0;

c_values = (40:20:140);
r_values = (40:20:140);


methods = {'deterministic','subspace_expected','uniform_sampling','near_optimal'};
out = run_different_number_of_c_and_r(in,methods,c_values,r_values);


save('./output/compare_abalone20_sigma_5')

c_values_plot;
saveas(gcf,'./plots/c_plots_abalone_20_sigma_5','fig');
export_fig(gcf,'./plots/c_plots_abalone_20_sigma_5.pdf');
close all;