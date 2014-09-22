load 'protein_data';
in.A = normalize_kernel_data(X);
clear X;

in.k = 20;
in.c = 50;
in.r = 50;
in.q = 10;

methods = {'subspace_expected','deterministic'};
p_values = (10:49);

out = run_dataset_different_p(in,methods,p_values);

save('./output/SNPs10')

p_values_plot_deterministic;
export_fig(gcf,'./plots/p_plots_SNPs10_deterministic.pdf');
close all;

p_values_plot_subspace_expected;
export_fig(gcf,'./plots/p_plots_SNPs10_subspace_expected.pdf');
close all;