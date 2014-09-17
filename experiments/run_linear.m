%run cran_field k = 5

load 'cranfield-term-doc'
in.A = normalize_kernel_data(A_cran);
in.c_and_r = [50,60,70,80,90,100;50,60,70,80,90,100];
in.k = 5;
in.q = 2;
in.p = [k+5,k+10,k+15,k+20];
in.methods = {'uniform_sampling','subspace_expeceted'};
fpinrtf('run cranfield-term-doc dataset with k=%d\n',in.k);
output_cranfield_5 = rundataset(in);