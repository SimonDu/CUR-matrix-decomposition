in.A = read_snap_data('CA-GrQc.txt');

in.k = 10;
in.p = 10;
in.q = 10;

in.sigma_k = 1;
in.froerr = 1;
in.froerr_k = 1;
in.specerr = 0;
in.specerr_k = 0;

c_values = (20:20:140);
r_values = 2*c_values;

in.adaptive = 1;

methods = {'deterministic','randomized_unweighted','subspace_expected','subspace_approxlevscore_gaussian','uniform_sampling','near_optimal'};
out = run_different_number_of_c_and_r(in,methods,c_values,r_values);


save('./output/compare_GR10_adaptive')

c_values_plot;
saveas(gcf,'./plots/c_plots_GR10_adaptive','fig');
export_fig(gcf,'./plots/c_plots_GR10_adaptive.pdf');
close all;