load 'snpsdata'
A = A(1:(end-3), 1:size(A,2)-3); % The last three columns and rows are all 0s
in.A = normalize_kernel_data(A);
in.A = in.A * in.A';
clear A;

in.k = 5;
in.c = 46;
in.r = 46;
in.q = 15;

methods = {'subspace_expected','deterministic'};
p_values = (5:20);

s = svds(in.A,p_values(end));
plot(p_values,s(p_values(1):end)./s(p_values(1)));
title('singular value decay');
export_fig('./plots/decay_SNPs5.pdf');
close all;


out = run_dataset_different_p(in,methods,p_values);

save('./output/SNPs5')

p_values_plot_deterministic;
export_fig(gcf,'./plots/p_plots_SNPs5_deterministic.pdf');
close all;

p_values_plot_subspace_expected;
export_fig(gcf,'./plots/p_plots_SNPs5_subspace_expected.pdf');
close all;