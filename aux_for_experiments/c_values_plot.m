A = in.A;
k = in.k;
[U,S,V] = svds(A,k);
fro_A_A_k = norm(A-U*S*V','fro');
spec_A_A_k = svds(A-U*S*V',1);



figure('position',[.1 .1 1300 200])
%sigma_k_plot
if(in.sigma_k)
    subplot(1,5,1);
    title('kth Singular Value Ratio','FontSize',15);
    hold on;
    
%     plot(c_values,out.randomized_unweighted.sigma_k./S(k,k),'Color','y','LineStyle','--');
    plot(c_values,out.subspace_expected.sigma_k./S(k,k),'-bo','LineWidth',2);
%     plot(c_values,out.subspace_approxlevscore_gaussian.sigma_k./S(k,k),'Color','b','LineStyle','--');
    plot(c_values,out.uniform_sampling.sigma_k./S(k,k),'-g+','LineWidth',2);
    plot(c_values,out.near_optimal.sigma_k./S(k,k),'-ks','LineWidth',2);
    plot(c_values,out.deterministic.sigma_k./S(k,k),'-rx','LineWidth',2);
    h_legend = legend('RandLeverage','RandUniform','NearOptimal','DetUCS','Location','southeast');
    set(h_legend,'FontSize',10);
    xlabel('value of c','FontSize',15);
    ylabel('\sigma_k(CUR)/\sigma_k(A)','FontSize',15);
    xlim([c_values(1),c_values(end)]);
    set(gca,'FontSize',12);
end

%froerr-plot
if(in.froerr)
    subplot(1,5,2);
    title('Frobenius Norm Error','FontSize',15);
    hold on;
    
%     plot(c_values,out.randomized_unweighted.froerr./fro_A_A_k,'Color','r','LineStyle','--');
    plot(c_values,out.subspace_expected.froerr./fro_A_A_k,'-bo','LineWidth',2);
%     plot(c_values,out.subspace_approxlevscore_gaussian.froerr./fro_A_A_k,'Color','b','LineStyle','--');
    plot(c_values,out.uniform_sampling.froerr./fro_A_A_k,'-g+','LineWidth',2);
    plot(c_values,out.near_optimal.froerr./fro_A_A_k,'-ks','LineWidth',2);
    plot(c_values,out.deterministic.froerr./fro_A_A_k,'-rx','LineWidth',2);
    h_legend = legend('RandLeverage','RandUniform','NearOptimal','DetUCS','Location','northeast');
    set(h_legend,'FontSize',10);
    xlabel('value of c','FontSize',15);
    ylabel('||A-CUR||_F/||A-A_k||_F','FontSize',15);
    xlim([c_values(1),c_values(end)]);
    set(gca,'FontSize',12);
end

%froerr_k
if(in.froerr_k)
    subplot(1,5,3);
    title('Rank k F-Norm Error','FontSize',15);
    hold on;
    
%     plot(c_values,out.randomized_unweighted.froerr_k./fro_A_A_k,'Color','r','LineStyle','--');
    plot(c_values,out.subspace_expected.froerr_k./fro_A_A_k,'-bo','LineWidth',2);
%     plot(c_values,out.subspace_approxlevscore_gaussian.froerr_k./fro_A_A_k,'Color','b','LineStyle','--');
    plot(c_values,out.uniform_sampling.froerr_k./fro_A_A_k,'-g+','LineWidth',2);
    plot(c_values,out.near_optimal.froerr_k./fro_A_A_k,'-ks','LineWidth',2);
    plot(c_values,out.deterministic.froerr_k./fro_A_A_k,'-rx','LineWidth',2);
    h_legend = legend('RandLeverage','RandUniform','NearOptimal','DetUCS','Location','northeast');
    set(h_legend,'FontSize',10);
    xlabel('value of c','FontSize',15);
    ylabel('||A-CUR_k||_F/||A-A_k||_F','FontSize',15);
    xlim([c_values(1),c_values(end)]);
    set(gca,'FontSize',12);
end

%specerr
if(in.specerr)
    subplot(1,5,4);
    title('Spectral Norm Error','FontSize',15);
    hold on;
    
%     plot(c_values,out.randomized_unweighted.specerr./spec_A_A_k,'Color','r','LineStyle','--');
    plot(c_values,out.subspace_expected.specerr./spec_A_A_k,'-bo','LineWidth',2);
%     plot(c_values,out.subspace_approxlevscore_gaussian.specerr./spec_A_A_k,'Color','b','LineStyle','--');
    plot(c_values,out.uniform_sampling.specerr./spec_A_A_k,'-g+','LineWidth',2);
    plot(c_values,out.near_optimal.specerr./spec_A_A_k,'-ks','LineWidth',2);
    plot(c_values,out.deterministic.specerr./spec_A_A_k,'-rx','LineWidth',2);
    h_legend = legend('RandLeverage','RandUniform','NearOptimal','DetUCS','Location','northeast');
    set(h_legend,'FontSize',10);
    xlabel('value of c','FontSize',15);
    ylabel('||A-CUR||_2/||A-A_k||_2','FontSize',15);
    xlim([c_values(1),c_values(end)]);
    set(gca,'FontSize',12);
end

%speccerr_k
if(in.specerr_k)
    subplot(1,5,5);
    title('Rank k Spectral Norm Error','FontSize',15);
    hold on;
    
%     plot(c_values,out.randomized_unweighted.specerr_k./spec_A_A_k,'Color','r','LineStyle','--');
    plot(c_values,out.subspace_expected.specerr_k./spec_A_A_k,'-bo','LineWidth',2);
%     plot(c_values,out.subspace_approxlevscore_gaussian.specerr_k./spec_A_A_k,'Color','b','LineStyle','--');
    plot(c_values,out.uniform_sampling.specerr_k./spec_A_A_k,'-g+','LineWidth',2);
    plot(c_values,out.near_optimal.specerr_k./spec_A_A_k,'-ks','LineWidth',2);
    plot(c_values,out.deterministic.specerr_k./spec_A_A_k,'-rx','LineWidth',2);
    h_legend = legend('RandLeverage','RandUniform','NearOptimal','DetUCS','Location','northeast');
    set(h_legend,'FontSize',10);
    xlabel('value of c','FontSize',15);
    ylabel('||A-CUR_k||_2/||A-A_k||_2','FontSize',15);
    xlim([c_values(1),c_values(end)]);
    set(gca,'FontSize',12);
end