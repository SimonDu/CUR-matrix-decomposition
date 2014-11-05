% Reuters data with target rank = 20
%load data
load('Reuters21578.mat');
in.A = fea;

in.k = 20;
in.p = 20;
in.q = 10;

in.sigma_k = 1;
in.froerr = 1;
in.froerr_k = 1;
in.specerr = 1;
in.specerr_k = 1;

in.adaptive = 0;

c_values = (40:20:140);
r_values = (40:20:140);


methods = {'deterministic','subspace_expected','uniform_sampling','near_optimal'};
out = run_different_number_of_c_and_r(in,methods,c_values,r_values);


save('./output/compare_reuters10')

c_values_plot;
saveas(gcf,'./plots/c_plots_reuters10','fig');
export_fig(gcf,'./plots/c_plots_reuters10.pdf');
close all;