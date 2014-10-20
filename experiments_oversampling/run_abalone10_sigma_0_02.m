in.sigma = 0.05;
load abalone_dataset;
X = generate_distance_matrix(abaloneInputs');
in.A = generate_RBF_kernel(X, in.sigma);
clear X;

in.k = 10;
in.c = 50;
in.r = 50;
in.q = 20;

in.sigma_k = 1;
in.froerr = 1;
in.froerr_k = 1;
in.specerr = 1;
in.specerr_k = 1;

in.adaptive = 0;

methods = {'subspace_expected','deterministic'};
p_values = (10:20);

s = svds(in.A,p_values(end));
plot(p_values,s(p_values(1):end)./s(p_values(1)));
title('singular value decay');
export_fig('./plots/decay_abalone10_sigma_0_02.pdf');
close all;

out = run_dataset_different_p(in,methods,p_values);

save('./output/p_plots_abalone10_sigma_0_02')

p_values_plot_subspace_expected;
saveas(gcf,'./plots/p_plots_abalone10_sigma_0_02_subspace_expected','fig');
export_fig(gcf,'./plots/p_plots_abalone10_sigma_0_02_subspace_expected.pdf');
close all;

p_values_plot_deterministic;
saveas(gcf,'./plots/p_plots_abalone10_sigma_0_02_deterministic','fig');
export_fig(gcf,'./plots/p_plots_abalone10_sigma_0_02_deterministic.pdf');
close all;
