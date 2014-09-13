load('/Users/simonshaoleidu/Programming/github/CX and CUR matrix decomposition/datasets/cranfield-term-doc.mat');
in_cell = {};
out_cell = {};
k=5;
number_of_c_and_r = [50,50];


[U,S,V] = svds(A_cran,5);
s=svds(A_cran,5);
A_k = U*S*V';
froerr = norm(A_cran-A_k,'fro');
specerr = norm(A_cran-A_k);


for i = 5:30
    in_cell{i}.A = A_cran;
    in_cell{i}.k=k;
    in_cell{i}.p=i;
    in_cell{i}.q=1;
    in_cell{i}.number_of_c_and_r = number_of_c_and_r;
    out_cell{i} = deterministic(in_cell{i});
    disp(i);
end

hold;
sigma_k_cell=[];
for i=5:30
    sigma_k_cell(i-4) = out_cell{i}.sigma_k / s(5);
end
plot((5:30),sigma_k_cell);


hold;
fro_cell=[];
for i=5:30
    fro_cell(i-4) = out_cell{i}.froerr(1) /froerr;
end
plot((5:30),fro_cell);

hold;
fro_k_cell=[];
for i=5:30
    fro_k_cell(i-4) = out_cell{i}.froerr(2)/froerr;
end
plot((5:30),fro_k_cell);

hold;
fro_cell=[];
for i=5:30
    fro_cell(i-4) = out_cell{i}.specerr(1)/specerr;
end
plot((5:30),fro_cell);

hold;
fro_k_cell=[];
for i=5:30
    fro_k_cell(i-4) = out_cell{i}.specerr(2)/specerr;
end
plot((5:30),fro_k_cell);