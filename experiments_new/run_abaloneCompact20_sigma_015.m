cutoffmultiplier = 3;
in.sigma = 0.15; 
in.cutoff = cutoffmultiplier*in.sigma;
in.d = 8;
load abalone_dataset;
X = generate_distance_matrix(abaloneInputs');
in.A = generate_compact_RBF_kernel(X, in.sigma, in.d, in.cutoff);
clear X;

in.k = 20;
in.c = 150;
in.r = 150;
in.q = 15;


methods = {'subspace_expected','deterministic'};
p_values = (20:40);

s = svds(in.A,p_values(end));
plot(s(p_values(1):end)./s(p_values(1)));
title('singular value decay');
xlim([p_values(1),p_values(end)]);
export_fig('./plots/decay_abaloneCompact20_sigma_015_q_15.pdf');
close all;

out = run_dataset_different_p(in,methods,p_values);

save('./output/abaloneCompact20_sigma_015_q_15')

p_values_plot_deterministic;
export_fig(gcf,'./plots/p_plots_abaloneCompact20_sigma_015_deterministic_q_15.pdf');
close all;

p_values_plot_subspace_expected;
export_fig(gcf,'./plots/p_plots_abaloneCompact20_sigma_015_subspace_expected_q_15.pdf');
close all;