
A = in.A;
k = in.k;
[U,S,V] = svds(A,k);
fro_A_A_k = norm(A-U*S*V','fro');
spec_A_A_k = svds(A-U*S*V',1);

figure;
%sigma_k_plot
if(in.sigma_k)
    subplot(2,3,1);
    title('Kth Singular Value Ratio');
    hold on;
    plot(p_values,out.randomized_unweighted.sigma_k./S(k,k),'r');
    h_legend = legend('randomized unweighted','Location','northeast');
    set(h_legend,'FontSize',5);
    xlabel('value of p','FontSize',15);
    ylabel('\sigma_k(CUR)/\sigma_k(A)','FontSize',15);
    xlim([p_values(1),p_values(end)]);
end

%froerr-plot
if(in.froerr)
    subplot(2,3,2);
    title('Frobenius Norm Error');
    hold on;
    plot(p_values,out.randomized_unweighted.froerr./fro_A_A_k,'r');
    h_legend = legend('randomized unweighted','Location','northeast');
    set(h_legend,'FontSize',5);
    xlabel('value of p','FontSize',15);
    ylabel('||A-CUR||_F/||A-A_k||_F','FontSize',15);
    xlim([p_values(1),p_values(end)]);
end

%froerr_k
if(in.froerr_k)
    subplot(2,3,3);
    title('Rank K Frobenius Norm Error');
    hold on;
    plot(p_values,out.randomized_unweighted.froerr_k./fro_A_A_k,'r');
    h_legend = legend('randomized unweighted','Location','northeast');
    set(h_legend,'FontSize',5);
    xlabel('value of p','FontSize',15);
    ylabel('||A-CUR_k||_F/||A-A_k||_F','FontSize',15);
    xlim([p_values(1),p_values(end)]);
end

%specerr
if(in.specerr)
    subplot(2,3,4);
    title('Spectral Norm Error');
    hold on;
    plot(p_values,out.randomized_unweighted.specerr./spec_A_A_k,'r');
    h_legend = legend('randomized unweighted','Location','northeast');
    set(h_legend,'FontSize',5);
    xlabel('value of p','FontSize',15);
    ylabel('||A-CUR||_2/||A-A_k||_2','FontSize',15);
    xlim([p_values(1),p_values(end)]);
end

%speccerr_k
if(in.specerr_k)
    subplot(2,3,5);
    title('Rank K Spectral Norm Error');
    hold on;
    plot(p_values,out.randomized_unweighted.specerr_k./spec_A_A_k,'r');
    h_legend = legend('randomized unweighted','Location','northeast');
    set(h_legend,'FontSize',5);
    xlabel('value of c','FontSize',15);
    ylabel('||A-CUR_k||_2/||A-A_k||_2','FontSize',15);
    xlim([p_values(1),p_values(end)]);
end