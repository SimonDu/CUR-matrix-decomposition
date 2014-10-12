function A = read_dexter

B = load('./datasets/dextertestdata.sparse');
A = sparse(B(:,1),B(:,2),B(:,3));