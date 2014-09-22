load 'snpsdata'
A = A(1:(end-3), 1:size(A,2)-3); % The last three columns and rows are all 0s
in.A = normalize_kernel_data(A*A');
clear A;

in.k = 5;
in.c = floor(6*in.k*log(in.k));
in.r = floor(6*in.k*log(in.k));
in.q = 10;

methods = {'subspace_expected','deterministic'};
p_values = (5:floor(6*in.k*log(in.k))-1);

out = run_dataset_different_p(in,methods,p_values);

save('./output/SNPs5')

p_values_plot_deterministic;
export_fig(gcf,'./plots/p_plots_SNPs5_deterministic.pdf');
close all;

p_values_plot_subspace_expected;
export_fig(gcf,'./plots/p_plots_SNPs5_subspace_expected.pdf');
close all;