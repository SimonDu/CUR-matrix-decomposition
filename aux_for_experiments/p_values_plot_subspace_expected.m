
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
    plot(p_values,out.subspace_expected.sigma_k./S(k,k),'-bo','LineWidth',2);
    h_legend = legend('subspace');
    set(h_legend,'FontSize',12);
    xlabel('value of p','FontSize',15);
    ylabel('\sigma_k(CUR)/\sigma_k(A)','FontSize',15);
    xlim([p_values(1),p_values(end)]);
end

%froerr-plot
if(in.froerr)
    subplot(1,5,2);
    title('Frobenius Norm Error','FontSize',15);
    hold on;
    plot(p_values,out.subspace_expected.froerr./fro_A_A_k,'-bo','LineWidth',2);
    h_legend = legend('subspace');
    set(h_legend,'FontSize',12);
    xlabel('value of p','FontSize',15);
    ylabel('||A-CUR||_F/||A-A_k||_F','FontSize',15);
    xlim([p_values(1),p_values(end)]);
end

%froerr_k
if(in.froerr_k)
    subplot(1,5,3);
    title('Rank-k Frobenius Norm Error','FontSize',15);
    hold on;
    plot(p_values,out.subspace_expected.froerr_k./fro_A_A_k,'-bo','LineWidth',2);
    h_legend = legend('subspace');
    set(h_legend,'FontSize',12);
    xlabel('value of p','FontSize',15);
    ylabel('||A-CUR_k||_F/||A-A_k||_F','FontSize',15);
    xlim([p_values(1),p_values(end)]);
end

%specerr
if(in.specerr)
    subplot(1,5,4);
    title('Spectral Norm Error','FontSize',15);
    hold on;
    plot(p_values,out.subspace_expected.specerr./spec_A_A_k,'-bo','LineWidth',2);
    h_legend = legend('subspace');
    set(h_legend,'FontSize',12);
    xlabel('value of p','FontSize',15);
    ylabel('||A-CUR||_2/||A-A_k||_2','FontSize',15);
    xlim([p_values(1),p_values(end)]);
end

%speccerr_k
if(in.specerr_k)
    subplot(1,5,5);
    title('Rank-k Spectral Norm Error','FontSize',15);
    hold on;
    plot(p_values,out.subspace_expected.specerr_k./spec_A_A_k,'-bo','LineWidth',2);
    h_legend = legend('subspace');
    set(h_legend,'FontSize',12);
    xlabel('value of c','FontSize',15);
    ylabel('||A-CUR_k||_2/||A-A_k||_2','FontSize',15);
    xlim([p_values(1),p_values(end)]);
end