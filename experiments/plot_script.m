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
        subspace_sigma_k=[subspace_sigma_k,subspace_output{i,j}.sigma_k];
    end
    subplot(2,3,1);
    plot(deterministic_sigma_k,'r');
    plot(subspace_sigma_k,'b');
    xlabel('number of columns and rows');
    xlim(c(1),c(end));
    %froerr_plot
%     deterministic_froerr=[];
%     subspace_froerr=[];
%     
%     subplot(2,3,2);
%     %froerr_k_plot
%     subplot(2,3,3);
%     %specerr_plot
%     subplot(2,3,4);
%     %specerr_k_plot
%     subplot(2,3,5)
end


% figure(1);
% for i=1:length(c)
%     sigma_k=[];
%     for j=1:length(p)
%         sigma_k=[sigma_k,deterministic_output{i,j}.sigma_k];
%     end
%     sigma_k = sigma_k./target_singular_value;
%     subplot(2,3,i);
%     plot(p,sigma_k,'r');
%     xlim([p(1) p(end)]);
% end
% 
% figure(2);
% for i=1:length(c)
%     froerrs=[];
%     for j=1:length(p)
%         froerrs=[froerrs,deterministic_output{i,j}.froerr(1)];
%     end
%     froerrs = froerrs./fro_A_A_k;
%     subplot(2,3,i);
%     plot(p,froerrs,'r');
%     xlim([p(1) p(end)]);
% end
% 
% figure(3);
% for i=1:length(c)
%     froerrs_k=[];
%     for j=1:length(p)
%         froerrs_k=[froerrs_k,deterministic_output{i,j}.froerr(2)];
%     end
%     froerrs_k = froerrs_k./fro_A_A_k;
%     subplot(2,3,i);
%     plot(p,froerrs_k,'r');
%     xlim([p(1) p(end)]);
% end
% 
% figure(4);
% for i=1:length(c)
%     specerrs=[];
%     for j=1:length(p)
%         specerrs=[specerrs,deterministic_output{i,j}.specerr(1)];
%     end
%     specerrs = specerrs./spec_A_A_k;
%     subplot(2,3,i);
%     plot(p,specerrs,'r');
%     xlim([p(1) p(end)]);
% end
% 
% figure(5);
% for i=1:length(c)
%     specerrs_k=[];
%     for j=1:length(p)
%         specerrs_k=[specerrs_k,deterministic_output{i,j}.specerr(2)];
%     end
%     specerrs_k = specerrs_k./spec_A_A_k;
%     subplot(2,3,i);
%     plot(p,specerrs_k),'r';
%     xlim([p(1) p(end)]);
% end



