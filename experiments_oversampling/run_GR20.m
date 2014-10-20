in.A = read_snap_data('CA-GrQc.txt');

in.k = 20;
in.c = 100;
in.r = 100;
in.q = 10;

in.sigma_k = 1;
in.froerr = 1;
in.froerr_k = 1;
in.specerr = 1;
in.specerr_k = 1;

in.adaptive = 0;

methods = {'deterministic','subspace_expected'};
p_values = (20:30);


s = svds(in.A,p_values(end));
plot(p_values,s(p_values(1):end)./s(p_values(1)));
title('singular value decay');
export_fig('./plots/decay_GR20.pdf');
close all;

out = run_dataset_different_p(in,methods,p_values);

save('./output/p_plots_GR20')

p_values_plot_deterministic;
export_fig(gcf,'./plots/p_plots_GR20_deterministic.pdf');
close all;

p_values_plot_subspace_expected;
export_fig(gcf,'./plots/p_plots_GR20_subspace_expected.pdf');
close all;