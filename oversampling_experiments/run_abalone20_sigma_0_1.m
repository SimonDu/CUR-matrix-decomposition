%oversampling experiments for RBF kernel of abalone dataset with sigma = 0.1,
%target rank = 20
in.sigma = 0.1;
load abalone_dataset;
X = generate_distance_matrix(abaloneInputs');
in.A = generate_RBF_kernel(X, in.sigma);
clear X;

in.k = 20;
in.c = 80;
in.r = 80;
in.q = 5;

in.sigma_k = 1;
in.froerr = 1;
in.froerr_k = 1;
in.specerr = 1;
in.specerr_k = 1;

in.adaptive = 0;

methods = {'subspace_expected','deterministic'};
p_values = (20:40);

s = svds(in.A,p_values(end));
plot(p_values,s(p_values(1):end)./s(p_values(1)));
title('singular value decay');
export_fig('./plots/decay_abalone20_sigma_0_1.pdf');
close all;

out = run_dataset_different_p(in,methods,p_values);

save('./output/p_plots_abalone20_sigma_0_1')
out_01 = out;
% p_values_plot_deterministic;
% saveas(gcf,'./plots/p_plots_abalone20_sigma_0_1_deterministic','fig');
% export_fig(gcf,'./plots/p_plots_abalone20_sigma_0_1_deterministic.pdf');
% close all;
% 
% p_values_plot_subspace_expected;
% saveas(gcf,'./plots/p_plots_abalone20_sigma_0_1_subspace_expected','fig');
% export_fig(gcf,'./plots/p_plots_abalone20_sigma_0_1_subspace_expected.pdf');
% close all;

