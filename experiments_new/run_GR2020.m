in.A = read_snap_data('CA-GrQc.txt');

in.k = 20;
in.c = 80;
in.r = 80;
in.q = 10;

methods = {'subspace_expected','deterministic'};
p_values = (20:79);

out = run_dataset_different_p(in,methods,p_values);

save('./output/GR20')

p_values_plot_deterministic;
export_fig(gcf,'./plots/p_plots_GR20_deterministic.pdf');
close all;

p_values_plot_subspace_expected;
export_fig(gcf,'./plots/p_plots_GR20_subspace_expected.pdf');
close all;