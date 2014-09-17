%set up
k=in.k;
A = in.A;
[U,S,V] = svds(A,k);
target_singular_value = S(k,k);
fro_A_A_k = norm(A-U*S*V','fro');
spec_A_A_k = svds(A-U*S*V',1);
%p,c for plot
p = (k:k+20);
c = (50:10:100);

%optimals for determinisitc algorithms
optimal_deterministic = {};
for i=1:length(c)
    optimal_deterministic{i} = find_optimal_p(deterministic_output,i);
end
%optimals for subspace algorithms
optimal_subspace = {};
for i=1:length(c)
    optimal_subspace{i} = find_optimal_p(subspace_output,i);
end

%plot sigma_k for different algorithms
sigma_k_deterministic=[];
for j=1:length(c)
    sigma_k_deterministic = [sigma_k_deterministic,optimal_deterministic{j}.sigma_k];
end
sigma_k_deterministic = sigma_k_deterministic./target_singular_value;
sigma_k_subspace=[];
for j=1:length(c)
    sigma_k_subspace = [sigma_k_subspace,optimal_subspace{j}.sigma_k];
end
sigma_k_subspace = sigma_k_subspace./target_singular_value;
sigma_k_uniform=[];
for j=1:length(c)
    sigma_k_uniform = [sigma_k_uniform,uniform_output{j,end}.sigma_k];
end
sigma_k_uniform = sigma_k_uniform./target_singular_value;
subplot(2,3,1);
hold on;
plot(c,sigma_k_deterministic,'r');
plot(c,sigma_k_subspace,'b');
plot(c,sigma_k_uniform,'g');
legend('deterministic','subspace','uniform');
xlabel('number of columns and rows');
title('ratio of kth singular value of CUR and that of A');
xlim([c(1) c(end)]);

%plot frobenius norm error for different algorithms
froerr_deterministic=[];
for j=1:length(c)
    froerr_deterministic = [froerr_deterministic,optimal_deterministic{j}.froerr];
end
froerr_deterministic = froerr_deterministic./fro_A_A_k;
froerr_subspace=[];
for j=1:length(c)
    froerr_subspace = [froerr_subspace,optimal_subspace{j}.froerr];
end
froerr_subspace = froerr_subspace./fro_A_A_k;
froerr_uniform=[];
for j=1:length(c)
    froerr_uniform = [froerr_uniform,mean(uniform_output{j,end}.froerr(1,:))];
end
froerr_uniform = froerr_uniform./fro_A_A_k;
subplot(2,3,2);
hold on;
plot(c,froerr_deterministic,'r');
plot(c,froerr_subspace,'b');
plot(c,froerr_uniform,'g');
legend('deterministic','subspace','uniform');
xlabel('number of columns and rows');
title('relative frobenius norm reconstruction error');
xlim([c(1) c(end)]);

%plot truncated-k frobenius norm error for different algorithms
froerr_k_deterministic=[];
for j=1:length(c)
    froerr_k_deterministic = [froerr_k_deterministic,optimal_deterministic{j}.froerr_k];
end
froerr_k_deterministic = froerr_k_deterministic./fro_A_A_k;
froerr_k_subspace=[];
for j=1:length(c)
    froerr_k_subspace = [froerr_k_subspace,optimal_subspace{j}.froerr_k];
end
froerr_k_subspace = froerr_k_subspace./fro_A_A_k;
froerr_k_uniform=[];
for j=1:length(c)
    froerr_k_uniform = [froerr_k_uniform,mean(uniform_output{j,end}.froerr(2,:))];
end
froerr_k_uniform = froerr_k_uniform./fro_A_A_k;
subplot(2,3,3);
hold on;
plot(c,froerr_k_deterministic,'r');
plot(c,froerr_k_subspace,'b');
plot(c,froerr_k_uniform,'g');
legend('deterministic','subspace','uniform');
xlabel('number of columns and rows');
title('relative frobenius norm rank-k reconstruction error');
xlim([c(1) c(end)]);

%plot spectral norm error for different algorithms
specerr_deterministic=[];
for j=1:length(c)
    specerr_deterministic = [specerr_deterministic,optimal_deterministic{j}.specerr];
end
specerr_deterministic = specerr_deterministic./spec_A_A_k;
specerr_subspace=[];
for j=1:length(c)
    specerr_subspace = [specerr_subspace,optimal_subspace{j}.specerr];
end
specerr_subspace = specerr_subspace./spec_A_A_k;
specerr_uniform=[];
for j=1:length(c)
    specerr_uniform = [specerr_uniform,mean(uniform_output{j,end}.specerr(1,:))];
end
specerr_uniform = specerr_uniform./spec_A_A_k;
subplot(2,3,4);
hold on;
plot(c,specerr_deterministic,'r');
plot(c,specerr_subspace,'b');
plot(c,specerr_uniform,'g');
legend('deterministic','subspace','uniform');
xlabel('number of columns and rows');
title('relative spectral norm reconstruction error');
xlim([c(1) c(end)]);

%plot truncated-k spectral norm error for different algorithms
specerr_k_deterministic=[];
for j=1:length(c)
    specerr_k_deterministic = [specerr_k_deterministic,optimal_deterministic{j}.specerr_k];
end
specerr_k_deterministic = specerr_k_deterministic./spec_A_A_k;
specerr_k_subspace=[];
for j=1:length(c)
    specerr_k_subspace = [specerr_k_subspace,optimal_subspace{j}.specerr_k];
end
specerr_k_subspace = specerr_k_subspace./spec_A_A_k;
specerr_k_uniform=[];
for j=1:length(c)
    specerr_k_uniform = [specerr_k_uniform,mean(uniform_output{j,end}.specerr(2,:))];
end
specerr_k_uniform = specerr_k_uniform./spec_A_A_k;
subplot(2,3,5);
hold on;
plot(c,specerr_k_deterministic,'r');
plot(c,specerr_k_subspace,'b');
plot(c,specerr_k_uniform,'g');
legend('deterministic','subspace','uniform');
xlabel('number of columns and rows');
title('relative spectral norm rank-k reconstruction error');
xlim([c(1) c(end)]);

