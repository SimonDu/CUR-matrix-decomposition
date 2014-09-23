function A = readdexter

B = load('dataset/dexter_train-edit');
A = sparse(B(:,1),B(:,2),B(:,3));