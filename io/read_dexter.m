function A = read_dexter

B = load('../datasets/dexter_train-edit');
A = sparse(B(:,1),B(:,2),B(:,3));