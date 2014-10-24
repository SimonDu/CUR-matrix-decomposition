s1 = 1:100;
s2 = 1:100;
for i = 2:100
    s1(i) = 0.965*s1(i-1);
    s2(i) = 0.996*s2(i-1);
end



hold on;
plot(s1,'-b','LineWidth',2)
plot(s2,'-r','LineWidth',2)
xlabel('i','FontSize',20)
ylabel('\sigma_i','FontSize',20)
h_legend = legend('Fast Singular Value Decay','Slow Singular Value Decay');
set(h_legend,'FontSize',15,'Location','northeast');
set(gca,'FontSize',20);

saveas(gcf,'sigmas','fig');
export_fig(gcf,'sigmas.pdf');
close all;