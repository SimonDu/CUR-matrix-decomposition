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

figure(6);
for i=1:length(c)
    sigma_k=[];
    for j=1:length(p)
        sigma_k=[sigma_k,subspace_output{i,j}.sigma_k];
    end
    sigma_k = sigma_k./target_singular_value;
    subplot(2,3,i);
    plot(p,sigma_k);
    xlim([p(1) p(end)]);
end

figure(6);
for i=1:length(c)
    froerrs=[];
    for j=1:length(p)
        froerrs=[froerrs,subspace_output{i,j}.froerr(1)];
    end
    froerrs = froerrs./fro_A_A_k;
    subplot(2,3,i);
    plot(p,froerrs);
    xlim([p(1) p(end)]);
end

figure(8);
for i=1:length(c)
    froerrs_k=[];
    for j=1:length(p)
        froerrs_k=[froerrs_k,subspace_output{i,j}.froerr(2)];
    end
    froerrs_k = froerrs_k./fro_A_A_k;
    subplot(2,3,i);
    plot(p,froerrs_k);
    xlim([p(1) p(end)]);
end

figure(9);
for i=1:length(c)
    specerrs=[];
    for j=1:length(p)
        specerrs=[specerrs,subspace_output{i,j}.specerr(1)];
    end
    specerrs = specerrs./spec_A_A_k;
    subplot(2,3,i);
    plot(p,specerrs);
    xlim([p(1) p(end)]);
end

figure(10);
for i=1:length(c)
    specerrs_k=[];
    for j=1:length(p)
        specerrs_k=[specerrs_k,subspace_output{i,j}.specerr(2)];
    end
    specerrs_k = specerrs_k./spec_A_A_k;
    subplot(2,3,i);
    plot(p,specerrs_k);
    xlim([p(1) p(end)]);
end



