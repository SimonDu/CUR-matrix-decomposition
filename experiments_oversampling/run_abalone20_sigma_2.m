in.sigma = 2;
load abalone_dataset;
X = generate_distance_matrix(abaloneInputs');
in.A = generate_RBF_kernel(X, in.sigma);
clear X;

in.k = 20;
in.c = 150;
in.r = 150;
in.q = 10;

in.sigma_k = 1;
in.froerr = 1;
in.froerr_k = 1;
in.specerr = 0;
in.specerr_k = 0;

in.adaptive = 0;

methods = {'deterministic','subspace_expected','subspace_approxlevscore_gaussian','randomized_unweighted'};
p_values = (20:40);

s = svds(in.A,p_values(end));
plot(p_values,s(p_values(1):end)./s(p_values(1)));
title('singular value decay');
export_fig('./plots/decay_abalone20_sigma_2.pdf');
close all;

out = run_dataset_different_p(in,methods,p_values);

save('./output/abalone20_sigma_2')

p_values_plot_deterministic;
export_fig(gcf,'./plots/p_plots_abalone20_sigma_2_deterministic.pdf');
close all;

p_values_plot_subspace_expected;
export_fig(gcf,'./plots/p_plots_abalone20_sigma_2_subspace_expected.pdf');
close all;

p_values_plot_subspace_approxlevscore_gaussian;
export_fig(gcf,'./plots/p_plots_abalone20_sigma_2_subspace_approxlevscore_gaussian.pdf');
close all;

p_values_plot_randomized_unweighted;
export_fig(gcf,'./plots/p_plots_abalone20_sigma_2_randomized_unweighted.pdf');
close all;