A = rand(100);
[U,S,V]=svd(A);
for i =1:100
    S(i,i) = 0.99^(i);
end
k=10;
p=20;
hold on;
plot(diag(S),'LineWidth',2)
xlabel('i','FontSize',20)
ylabel('\sigma_i','FontSize',20)