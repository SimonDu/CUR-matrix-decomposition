A = in.A;
k = in.k;
[U,S,V] = svds(A,k);
fro_A_A_k = norm(A-U*S*V','fro');
spec_A_A_k = svds(A-U*S*V',1);

figure;
%sigma_k_plot
if(in.sigma_k)
    subplot(1,3,1);
    title('sigma_k');
    hold on;
    plot(c_values,out.deterministic.sigma_k./S(k,k),'r');
    plot(c_values,out.subspace_expected.sigma_k./S(k,k),'b');
    plot(c_values,out.subspace_approxlevscore_gaussian.sigma_k./S(k,k),'g');
    plot(c_values,out.randomized_unweighted.sigma_k./S(k,k),'y');
    plot(c_values,out.uniform_sampling.sigma_k./S(k,k),'--');
    legend('deterministic','exact subspace','gaussian subspace','randomized unweighted','uniform sampling');
    xlabel('value of c');
    xlim([c_values(1),c_values(end)]);
end

%froerr-plot
if(in.froerr)
    subplot(1,3,2);
    title('froerr');
    hold on;
    plot(c_values,out.deterministic.froerr./fro_A_A_k,'r');
    plot(c_values,out.subspace_expected.froerr./fro_A_A_k,'b');
    plot(c_values,out.subspace_approxlevscore_gaussian.froerr./fro_A_A_k,'g');
    plot(c_values,out.randomized_unweighted.froerr./fro_A_A_k,'y');
    plot(c_values,out.uniform_sampling.froerr./fro_A_A_k,'--');
    legend('deterministic','exact subspace','gaussian subspace','randomized unweighted','uniform sampling');
    xlabel('value of c');
    xlim([c_values(1),c_values(end)]);
end

%froerr_k
if(in.froerr_k)
    subplot(1,3,3);
    title('froerr-k');
    hold on;
    plot(c_values,out.deterministic.froerr_k./fro_A_A_k,'r');
    plot(c_values,out.subspace_expected.froerr_k./fro_A_A_k,'b');
    plot(c_values,out.subspace_approxlevscore_gaussian.froerr_k./fro_A_A_k,'g');
    plot(c_values,out.randomized_unweighted.froerr_k./fro_A_A_k,'y');
    plot(c_values,out.uniform_sampling.froerr_k./fro_A_A_k,'--');
    legend('deterministic','exact subspace','gaussian subspace','randomized unweighted','uniform_sampling');
    xlabel('value of c');
    xlim([c_values(1),c_values(end)]);
end

%specerr
if(in.specerr)
    subplot(2,3,4);
    title('specerr');
    hold on;
    plot(c_values,out.deterministic.specerr./spec_A_A_k,'r');
    plot(c_values,out.subspace_expected.specerr./spec_A_A_k,'b');
    plot(c_values,out.subspace_approxlevscore_gaussian.specerr./specerr_A_A_k,'g');
    plot(c_values,out.randomized_unweighted.specerr./specerr_A_A_k,'y');
    plot(c_values,out.uniform_sampling.specerr./spec_A_A_k,'--');
    legend('deterministic','exact subspace','gaussian subspace','randomized unweighted','uniform_sampling');
    xlabel('value of c');
    xlim([c_values(1),c_values(end)]);
end

%speccerr_k
if(in.specerr_k)
    subplot(2,3,5);
    title('specerr-k');
    hold on;
    plot(c_values,out.deterministic.specerr_k./spec_A_A_k,'r');
    plot(c_values,out.subspace_expected.specerr_k./spec_A_A_k,'b');
    plot(c_values,out.subspace_approxlevscore_gaussian.specerr+k./specerr_A_A_k,'g');
    plot(c_values,out.randomized_unweighted.specerr_k./specerr_A_A_k,'y');
    plot(c_values,out.uniform_sampling.specerr_k./spec_A_A_k,'--');
    legend('deterministic','exact subspace','gaussian subspace','randomized unweighted','uniform_sampling');
    xlabel('value of c');
    xlim([c_values(1),c_values(end)]);
end