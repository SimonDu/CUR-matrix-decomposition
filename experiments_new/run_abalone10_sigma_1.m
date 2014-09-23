in.sigma = 1;
load abalone_dataset;
X = generate_distance_matrix(abaloneInputs');
in.A = generate_RBF_kernel(X, in.sigma);
clear X;

in.k = 10;
in.c = 50;
in.r = 50;
in.q = 5;

methods = {'subspace_expected','deterministic'};
p_values = (10:49);

s = svds(in.A,p_values(end));
plot(s(p_values(1):end)./s(p_values(1)));
title('singular value decay');
export_fig('./plots/decay_abalone10_sigma_1.pdf');
close all;

out = run_dataset_different_p(in,methods,p_values);

save('./output/abalone10_sigma_1')

p_values_plot_deterministic;
export_fig(gcf,'./plots/p_plots_abalone10_sigma_1_deterministic.pdf');
close all;

p_values_plot_subspace_expected;
export_fig(gcf,'./plots/p_plots_abalone10_sigma_1_subspace_expected.pdf');
close all;