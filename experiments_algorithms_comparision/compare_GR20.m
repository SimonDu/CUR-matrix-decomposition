in.A = read_snap_data('CA-GrQc.txt');

in.k = 20;
in.p = 20;
in.q = 15;

in.sigma_k = 1;
in.froerr = 1;
in.froerr_k = 1;
in.specerr = 0;
in.specerr_k = 0;

c_values = (40:10:150);
r_values = (40:10:150);


methods = {'deterministic','subspace_expected','subspace_approxlevscore_gaussian','uniform_sampling'};
out = run_different_number_of_c_and_r(in,methods,c_values,r_values);


save('./output/compare_GR20')

c_values_plot;
export_fig(gcf,'./plots/c_plots_GR20.pdf');
close all;