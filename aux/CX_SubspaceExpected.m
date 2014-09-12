function [idx] = CX_SubspaceExpected(A,k,c)
%subspace sampling CUR

[~,n] = size(A);

p = zeros(n,1);

%truncated SVD
[~,~,Va]=svds(A,k);

%calculate the stat leverage score
for j = 1:n
    p(j) = 1/k*sum(Va(j,:).^2);
end


%choose min of c*p(j) and 1
prob = zeros(n,1);
for j = 1:n
    prob(j) = min(c*p(j),1);
end

%sample
idx = [];
pick = zeros(n,1);
%t = 1;
for j = 1:n
    pick(j) = randsample([0 1],1,true,[1-prob(j) prob(j)]);
    if pick(j) == 1
        %S(j,t) = 1;
        %D(t,t)  = 1/min(1,sqrt(c*p(j)));
        %t = t + 1;
        idx = [idx j];
    end
end