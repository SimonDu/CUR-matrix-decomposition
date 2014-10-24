load abalone_dataset;
X = generate_distance_matrix(abaloneInputs');
A5 = generate_RBF_kernel(X, 5);
A01 = generate_RBF_kernel(X, 0.1);
clear X;

k = 20;
p_values = (20:40);
[U5,S5,V5] = svds(A5,k);
[U01,S01,V01] = svds(A01,k);
fro_A_A_k_5 = norm(A5-U5*S5*V5','fro');
spec_A_A_k_5 = svds(A5-U5*S5*V5',1);
fro_A_A_k_01 = norm(A01-U01*S01*V01','fro');
spec_A_A_k_01 = svds(A01-U01*S01*V01',1);


%froerr
figure('position',[.1 .1 1300 200])
subplot(1,4,1);
title('\sigma = 5','FontSize',15);
hold on;
plot(p_values,sigma_5_deter.froerr./fro_A_A_k_5,'-rx','LineWidth',2);
h_legend = legend('DetUCS');
set(h_legend,'FontSize',12);
xlabel('value of p','FontSize',15);
ylabel('||A-CUR||_F/||A-A_k||_F','FontSize',15);
xlim([p_values(1),p_values(end)]);

subplot(1,4,2);
title('\sigma = 5','FontSize',15);
hold on;
plot(p_values,sigma_5_subspace.froerr./fro_A_A_k_5,'-bo','LineWidth',2);
h_legend = legend('RandLeverage');
set(h_legend,'FontSize',12);
xlabel('value of p','FontSize',15);
ylabel('||A-CUR||_F/||A-A_k||_F','FontSize',15);
xlim([p_values(1),p_values(end)]);


subplot(1,4,3);
title('\sigma = 0.1','FontSize',15);
hold on;
plot(p_values,sigma_01_deter.froerr./fro_A_A_k_01,'-rx','LineWidth',2);
h_legend = legend('DetUCS');
set(h_legend,'FontSize',12);
xlabel('value of p','FontSize',15);
ylabel('||A-CUR||_F/||A-A_k||_F','FontSize',15);
xlim([p_values(1),p_values(end)]);


subplot(1,4,4);
title('\sigma = 0.1','FontSize',15);
hold on;
plot(p_values,sigma_01_subspace.froerr./fro_A_A_k_01,'-bo','LineWidth',2);
h_legend = legend('RandLeverage');
set(h_legend,'FontSize',12);
xlabel('value of p','FontSize',15);
ylabel('||A-CUR||_F/||A-A_k||_F','FontSize',15);
xlim([p_values(1),p_values(end)]);

saveas(gcf,'./p_plots/froerr','fig');
export_fig(gcf,'./p_plots/froerr.pdf');
close all;


%specerr
figure('position',[.1 .1 1300 200])
subplot(1,4,1);
title('\sigma = 5','FontSize',15);
hold on;
plot(p_values,sigma_5_deter.specerr./spec_A_A_k_5,'-rx','LineWidth',2);
h_legend = legend('DetUCS');
set(h_legend,'FontSize',12);
xlabel('value of p','FontSize',15);
ylabel('||A-CUR||_2/||A-A_k||_2','FontSize',15);
xlim([p_values(1),p_values(end)]);

subplot(1,4,2);
title('\sigma = 5','FontSize',15);
hold on;
plot(p_values,sigma_5_subspace.specerr./spec_A_A_k_5,'-bo','LineWidth',2);
h_legend = legend('RandLeverage');
set(h_legend,'FontSize',12);
xlabel('value of p','FontSize',15);
ylabel('||A-CUR||_2/||A-A_k||_2','FontSize',15);
xlim([p_values(1),p_values(end)]);


subplot(1,4,3);
title('\sigma = 0.1','FontSize',15);
hold on;
plot(p_values,sigma_01_deter.specerr./spec_A_A_k_01,'-rx','LineWidth',2);
h_legend = legend('DetUCS');
set(h_legend,'FontSize',12);
xlabel('value of p','FontSize',15);
ylabel('||A-CUR||_2/||A-A_k||_2','FontSize',15);
xlim([p_values(1),p_values(end)]);


subplot(1,4,4);
title('\sigma = 0.1','FontSize',15);
hold on;
plot(p_values,sigma_01_subspace.specerr./spec_A_A_k_01,'-bo','LineWidth',2);
h_legend = legend('RandLeverage');
set(h_legend,'FontSize',12,'Location','southeast');
xlabel('value of p','FontSize',15);
ylabel('||A-CUR||_2/||A-A_k||_2','FontSize',15);
xlim([p_values(1),p_values(end)]);

saveas(gcf,'./p_plots/specerr','fig');
export_fig(gcf,'./p_plots/specerr.pdf');
close all;


%sigmak
figure('position',[.1 .1 1300 200])
subplot(1,4,1);
title('\sigma = 5','FontSize',15);
hold on;
plot(p_values,sigma_5_deter.sigma_k./S5(k,k),'-rx','LineWidth',2);
h_legend = legend('DetUCS');
set(h_legend,'FontSize',12);
xlabel('value of p','FontSize',15);
ylabel('\sigma_k(CUR)/\sigma_k(A)','FontSize',15);
xlim([p_values(1),p_values(end)]);

subplot(1,4,2);
title('\sigma = 5','FontSize',15);
hold on;
plot(p_values,sigma_5_subspace.sigma_k./S5(k,k),'-bo','LineWidth',2);
h_legend = legend('RandLeverage');
set(h_legend,'FontSize',12);
xlabel('value of p','FontSize',15);
ylabel('\sigma_k(CUR)/\sigma_k(A)','FontSize',15);
xlim([p_values(1),p_values(end)]);


subplot(1,4,3);
title('\sigma = 0.1','FontSize',15);
hold on;
plot(p_values,sigma_01_deter.sigma_k./S01(k,k),'-rx','LineWidth',2);
h_legend = legend('DetUCS');
set(h_legend,'FontSize',12);
xlabel('value of p','FontSize',15);
ylabel('\sigma_k(CUR)/\sigma_k(A)','FontSize',15);
xlim([p_values(1),p_values(end)]);


subplot(1,4,4);
title('\sigma = 0.1','FontSize',15);
hold on;
plot(p_values,sigma_01_subspace.sigma_k./S01(k,k),'-bo','LineWidth',2);
h_legend = legend('RandLeverage');
set(h_legend,'FontSize',12);
xlabel('value of p','FontSize',15);
ylabel('\sigma_k(CUR)/\sigma_k(A)','FontSize',15);
xlim([p_values(1),p_values(end)]);

saveas(gcf,'./p_plots/sigma_k','fig');
export_fig(gcf,'./p_plots/sigma_k.pdf');
close all;
