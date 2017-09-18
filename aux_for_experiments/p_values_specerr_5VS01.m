figure('position',[.1 .1 1300 200]);
subplot(1,4,1)
title('\sigma = 5','FontSize',15);
hold on;
plot(p_values,out_5.deterministic.specerr./fro_A_A_k,'-rx','LineWidth',2);
h_legend = legend('StableCUR');
set(h_legend,'FontSize',12);
xlabel('value of p','FontSize',15);
ylabel('||A-CUR||_2/||A-A_k||_2','FontSize',15);
xlim([p_values(1),p_values(end)]);

subplot(1,4,2)
title('\sigma = 5','FontSize',15);
hold on;
plot(p_values,out_5.subspace_expected.specerr./fro_A_A_k,'-bx','LineWidth',2);
h_legend = legend('RandLeverage');
set(h_legend,'FontSize',12);
xlabel('value of p','FontSize',15);
ylabel('||A-CUR||_2/||A-A_k||_2','FontSize',15);
xlim([p_values(1),p_values(end)]);

subplot(1,4,3)
title('\sigma = 0.1','FontSize',15);
hold on;
plot(p_values,out_01.deterministic.specerr./fro_A_A_k,'-rx','LineWidth',2);
h_legend = legend('StableCUR');
set(h_legend,'FontSize',12);
xlabel('value of p','FontSize',15);
ylabel('||A-CUR||_2/||A-A_k||_2','FontSize',15);
xlim([p_values(1),p_values(end)]);

subplot(1,4,4)
title('\sigma = 0.1','FontSize',15);
hold on;
plot(p_values,out_01.subspace_expected.specerr./fro_A_A_k,'-bx','LineWidth',2);
h_legend = legend('RandLeverage');
set(h_legend,'FontSize',12);
xlabel('value of p','FontSize',15);
ylabel('||A-CUR||_2/||A-A_k||_2','FontSize',15);
xlim([p_values(1),p_values(end)]);