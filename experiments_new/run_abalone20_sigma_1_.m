in.sigma = 1;
X = generate_abalone_dataset(in.sigma); 
in.A = generate_RBF_kernel(X, in.sigma);
clear X;

in.k = 20;
in.c = 80;
in.r = 80;
in.q = 30;

methods = {'subspace_expected','deterministic'};
p_values = (20:79);

out = run_dataset_different_p(in,methods,p_values);

save('./output/abalone20_sigma_1')

p_values_plot_deterministic;
export_fig(gcf,'./plots/p_plots_abalone20_sigma_1_deterministic.pdf');
close all;

p_values_plot_subspace_expected;
export_fig(gcf,'./plots/p_plots_abalone20_sigma_1_subspace_expected.pdf');
close all;