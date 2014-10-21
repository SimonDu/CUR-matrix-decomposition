in.sigma = 0.05;
load abalone_dataset;
X = generate_distance_matrix(abaloneInputs');
in.A = generate_RBF_kernel(X, in.sigma);
clear X;

in.k = 20;
in.c = 40;
in.r = 40;
in.q = 20;

in.sigma_k = 1;
in.froerr = 1;
in.froerr_k = 1;
in.specerr = 1;
in.specerr_k = 1;

in.adaptive = 0;

methods = {'subspace_expected','deterministic'};
p_values = (20:39);

s = svds(in.A,p_values(end));
plot(p_values,s(p_values(1):end)./s(p_values(1)));
title('Singular Value Decay','FontSize',15);
xlabel('value of p','FontSize',15);
ylabel('pth singular value','FontSize',15);
export_fig('./plots/decay_abalone10_sigma_2.pdf');
close all;

out = run_dataset_different_p(in,methods,p_values);

save('./output/p_plots_abalone10_sigma_01')

p_values_plot_subspace_expected;
saveas(gcf,'./plots/p_plots_abalone20_sigma_005_c_40_subspace_expected','fig');
export_fig(gcf,'./plots/p_plots_abalone20_sigma_005_c_40_subspace_expected.pdf');
close all;

p_values_plot_deterministic;
saveas(gcf,'./plots/p_plots_abalone20_sigma_005_c_40_deterministic','fig');
export_fig(gcf,'./plots/p_plots_abalone20_sigma_005_c_40_deterministic.pdf');
close all;
