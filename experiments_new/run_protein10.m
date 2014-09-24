load 'protein_data';
in.A = normalize_kernel_data(X);
clear A;

in.k = 10;
in.c = 110;
in.r = 110;
in.q = 15;

methods = {'subspace_expected','deterministic'};
p_values = (10:20);

s = svds(in.A,p_values(end));
plot(p_values,s(p_values(1):end)./s(p_values(1)));
title('singular value decay');
export_fig('./plots/decay_protein10.pdf');
close all;

out = run_dataset_different_p(in,methods,p_values);

save('./output/protein10')

p_values_plot_deterministic;
export_fig(gcf,'./plots/p_plots_protein10_deterministic.pdf');
close all;

p_values_plot_subspace_expected;
export_fig(gcf,'./plots/p_plots_protein10_subspace_expected.pdf');
close all;