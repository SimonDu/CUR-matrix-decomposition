function [idx] = CX_SubspaceExpected(Va,p,c)
%subspace sampling CUR

[n,~] = size(Va);

L = zeros(n,1);

%truncated SVD


%calculate the stat leverage score
for j = 1:n
    L(j) = 1/p*norm(Va(j,:))^2;
end


%choose min of c*p(j) and 1
prob = zeros(n,1);
for j = 1:n
    prob(j) = min(c*L(j),1);
end

%sample
idx = [];
pick = zeros(n,1);
%t = 1;
for j = 1:n
    pick(j) = randsample([0 1],1,true,[1-prob(j) prob(j)]);
    if pick(j) == 1
        idx = [idx j];
    end
end