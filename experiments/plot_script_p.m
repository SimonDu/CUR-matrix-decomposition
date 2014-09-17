%a script to discover see the effect of choosing different p

%set up
k=10;
A = in.A;
[U,S,V] = svds(A,k);
target_singular_value = S(k,k);
fro_A_A_k = norm(A-U*S*V','fro');
spec_A_A_k = svds(A-U*S*V',1);
%p,c for plot
p = (k:k+20);
c = (50:10:100);
%p for different determinisitc

for i=1:length(c) 
    figure(i);
    %sigma_k_plot
    deterministic_sigma_k=[];
    subspace_sigma_k=[];
    for j=1:length(p)
        deterministic_sigma_k=[deterministic_sigma_k,deterministic_output{i,j}.sigma_k];
        subspace_sigma_k=[subspace_sigma_k,mean(subspace_output{i,j}.sigma_k)];
    end
    subplot(2,3,1);
    hold on;
    plot(p,deterministic_sigma_k./target_singular_value,'r');
    plot(p,subspace_sigma_k./target_singular_value,'b');
    legend('deterministic','subspace');
    xlabel('number of p');
    xlim([p(1),p(end)]);
    
    %froerr_plot
    deterministic_froerr=[];
    subspace_froerr=[];
    for j=1:length(p)
        deterministic_froerr=[deterministic_froerr,deterministic_output{i,j}.froerr(1,:)];
        subspace_froerr=[subspace_froerr,mean(subspace_output{i,j}.froerr(1,:))];
    end
    subplot(2,3,2);
    hold on;
    plot(p,deterministic_froerr./fro_A_A_k,'r');
    plot(p,subspace_froerr./fro_A_A_k,'b');
    legend('deterministic','subspace');
    xlabel('number of p');
    xlim([p(1),p(end)]);
    
    %froerr_k_plot
    deterministic_froerr_k=[];
    subspace_froerr_k=[];
    for j=1:length(p)
        deterministic_froerr_k=[deterministic_froerr_k,deterministic_output{i,j}.froerr(2,:)];
        subspace_froerr_k=[subspace_froerr_k,mean(subspace_output{i,j}.froerr(2,:))];
    end
    subplot(2,3,3);
    hold on;
    plot(p,deterministic_froerr_k./fro_A_A_k,'r');
    plot(p,subspace_froerr_k./fro_A_A_k,'b');
    legend('deterministic','subspace');
    xlabel('number of p');
    xlim([p(1),p(end)]);
    
    %specerr_plot
    deterministic_specerr=[];
    subspace_specerr=[];
    for j=1:length(p)
        deterministic_specerr=[deterministic_specerr,deterministic_output{i,j}.specerr(1,:)];
        subspace_specerr=[subspace_specerr,mean(subspace_output{i,j}.specerr(1,:))];
    end
    subplot(2,3,4);
    hold on;
    plot(p,deterministic_specerr./spec_A_A_k,'r');
    plot(p,subspace_specerr./spec_A_A_k,'b');
    legend('deterministic','subspace');
    xlabel('number of p');
    xlim([p(1),p(end)]);
    
    %froerr_k_plot
    deterministic_specerr_k=[];
    subspace_specerr_k=[];
    for j=1:length(p)
        deterministic_specerr_k=[deterministic_specerr_k,deterministic_output{i,j}.specerr(2,:)];
        subspace_specerr_k=[subspace_specerr_k,mean(subspace_output{i,j}.specerr(2,:))];
    end
    subplot(2,3,5);
    hold on;
    plot(p,deterministic_specerr_k./spec_A_A_k,'r');
    plot(p,subspace_specerr_k./spec_A_A_k,'b');
    legend('deterministic','subspace');
    xlabel('number of p');
    xlim([p(1),p(end)]);
       
end



