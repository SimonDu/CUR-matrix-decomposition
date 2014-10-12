%load data
in.A = read_dexter;

in.k = 10;
in.p = 10;
in.q = 10;

in.sigma_k = 1;
in.froerr = 1;
in.froerr_k = 1;
in.specerr = 0;
in.specerr_k = 0;

c_values = (20:20:140);
r_values = (20:20:140);


methods = {'deterministic','randomized_unweighted','subspace_expected','subspace_approxlevscore_gaussian','uniform_sampling','near_optimal'};
out = run_different_number_of_c_and_r(in,methods,c_values,r_values);


save('./output/compare_dexter10')

c_values_plot;
saveas(gcf,'./plots/c_plots_dexter10','fig');
export_fig(gcf,'./plots/c_plots_dexter10.pdf');
close all;