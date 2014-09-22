load 'protein_data';
in.A = normalize_kernel_data(X);
clear A;

in.k = 10;
in.c = 50;
in.r = 50;
in.q = 10;

methods = {'subspace_expected','deterministic'};
p_values = (in.k:in.c-1);

out = run_dataset_different_p(in,methods,p_values);

save('./output/protein10')

p_values_plot_deterministic;
export_fig(gcf,'./plots/p_plots_protein10_deterministic.pdf');
close all;

p_values_plot_subspace_expected;
export_fig(gcf,'./plots/p_plots_protein10_subspace_expected.pdf');
close all;