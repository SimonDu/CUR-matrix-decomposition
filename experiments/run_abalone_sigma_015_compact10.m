% run Abalone
% Abalone, sigma = .15

in.sigma = .15; % when sigma = .15, k=20 captures < 18% of variance
in.cutoff = cutoffmultiplier*in.sigma;
in.d = 8;

load 'abalone_distance_matrix'; 
in.A = generate_compact_RBF_kernel(X, in.sigma, in.d, in.cutoff);
clear X;


in5.A = in.A;
in10.A = in.A;

%data description
k = 5;
in5.k = k;
in5.p = [k+5,k+10,k+15,k+20];
out5_description = data_description(in5);
disp('done description 5');
k = 10;
in10.k = k;
in10.p = [k+5,k+10,k+15,k+20];
out10_description = data_description(in10);
disp('done description 10');

%set parameters
c_and_r = [50,60,70,80,90,100;50,60,70,80,90,100];

in.k = 10;
in.q = 6;

%deterministic

deterministic_output = {};
subspace_output = {};
for i=1:length(c_and_r)
    in.number_of_c_and_r = c_and_r(:,i);
    for p=in.k:in.k+20      
        in.p = p;
        deterministic_output{i,p-in.k+1} = deterministic(in);
        fprintf('done c = %d, p = %d, time CUR = %d, time error = %d\n',...
            in.number_of_c_and_r(1),p,deterministic_output{i,p-in.k+1}.timings(1),...
            deterministic_output{i,p-in.k+1}.timings(2));    
    end
end

%subspace
subspace_output = {};
for i=1:length(c_and_r)
    in.number_of_c_and_r = c_and_r(:,i);
    for p=in.k:in.k+20      
        in.p = p;
        subspace_output{i,p-in.k+1} = subspace_expected(in);
        fprintf('done c = %d, p = %d, time CUR = %d, time error = %d\n',...
            in.number_of_c_and_r(1),p,subspace_output{i,p-in.k+1}.timings(1),...
            subspace_output{i,p-in.k+1}.timings(2));    
    end
end


%uniform
uniform_output = {};
for i=1:length(c_and_r)
    in.number_of_c_and_r = c_and_r(:,i);
    uniform_output{i,p-in.k+1} = uniform_sampling(in);
    fprintf('done c = %d, time CUR = %d, time error = %d\n',...
        in.number_of_c_and_r(1),uniform_output{i,p-in.k+1}.timings(1),...
        uniform_output{i,p-in.k+1}.timings(2));    
end

save('./output/abalone_sigma_015_compact10');
