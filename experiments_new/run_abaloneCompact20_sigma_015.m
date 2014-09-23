cutoffmultiplier = 3;
in.sigma = 0.15; 
in.cutoff = cutoffmultiplier*in.sigma;
in.d = 8;
load abalone_dataset;
X = generate_distance_matrix(abaloneInputs');
in.A = generate_compact_RBF_kernel(X, in.sigma, in.d, in.cutoff);
clear X;

in.k = 20;
in.c = 80;
in.r = 80;
in.q = 10;


methods = {'subspace_expected','deterministic'};
p_values = (20:79);

s = svds(in.A,p_values(end));
plot(s(p_values(1):end)./s(p_values(1)));
title('singular value decay');
xlim([p_values(1),p_values(end)]);
export_fig('./plots/decay_abaloneCompact20_sigma_015.pdf');
close all;

out = run_dataset_different_p(in,methods,p_values);

save('./output/abaloneCompact20_sigma_015')

p_values_plot_deterministic;
export_fig(gcf,'./plots/p_plots_abaloneCompact20_sigma_015_deterministic.pdf');
close all;

p_values_plot_subspace_expected;
export_fig(gcf,'./plots/p_plots_abaloneCompact20_sigma_015_subspace_expected.pdf');
close all;