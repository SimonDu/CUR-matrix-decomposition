% linear kernel for GR data set with target rank = 20

in.A = read_snap_data('CA-GrQc.txt');

in.k = 20;
in.p = 20;
in.q = 10;

in.sigma_k = 1;
in.froerr = 1;
in.froerr_k = 1;
in.specerr = 1;
in.specerr_k = 1;

c_values = (40:20:140);
r_values = (40:20:140);

in.adaptive = 0;

methods = {'deterministic','subspace_expected','uniform_sampling','near_optimal'};
out = run_different_number_of_c_and_r(in,methods,c_values,r_values);


save('./output/compare_GR20')

c_values_plot;
saveas(gcf,'./plots/c_plots_GR20','fig');
export_fig(gcf,'./plots/c_plots_GR20.pdf');
close all;