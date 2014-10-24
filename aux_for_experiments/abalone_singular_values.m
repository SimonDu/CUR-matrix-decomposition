load abalone_dataset;
X = generate_distance_matrix(abaloneInputs');
in.A = generate_RBF_kernel(X, 5);
s5=svds(in.A,50);
in.A = generate_RBF_kernel(X, 2);
s2=svds(in.A,50);
in.A = generate_RBF_kernel(X, 0.2);
s02=svds(in.A,50);
in.A = generate_RBF_kernel(X, 0.1);
s01=svds(in.A,50);

hold on
plot(20:40,s01(20:40)/s01(20),'-bo','LineWidth',2)
plot(20:40,s02(20:40)/s02(20),'-ro','LineWidth',2)
plot(20:40,s2(20:40)/s2(20),'-mo','LineWidth',2)
plot(20:40,s5(20:40)/s5(20),'-go','LineWidth',2)
xlim([20,40]);
xlabel('p','FontSize',20);
ylabel('\sigma_p(A)/\sigma_k(A)','FontSize',20);
h_legend = legend('\sigma = 0.1','\sigma = 0.2','\sigma = 2','\sigma = 5');
set(h_legend,'FontSize',20,'Location','southwest');
set(gca,'FontSize',20);

saveas(gcf,'abalone_singular_values','fig');
export_fig(gcf,'abalone_singular_values.pdf');
close all;