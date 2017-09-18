A = rand(100);
[U,S,V]=svd(A);
for i =1:100
    S(i,i) = 2^(-i);
end

in.A = U*S*V';
in.k = 10;
in.p = 1;
c = 2:100;
r = 2:100;


in.q = 1;

in.sigma_k = 0;
in.froerr = 1;
in.froerr_k = 0;
in.specerr = 0;
in.specerr_k = 0;
in.adaptive = 0;

iter = length(c);
for i = 1:iter
    in.c = c(i);
    in.r = r(i);
    out = deterministic(in);
    froerr(i) = out.froerr;

    [~,~,Va]=svds(in.A,in.p);
    cidx = MSelect(Va(:,1:in.p),in.p,c(i));
    C = in.A(:,cidx);

    [~,~,Va]=svds(in.A',in.p);
    ridx = MSelect(Va(:,1:in.p),in.p,c(i));
    R = in.A(ridx,:);
    U = (C\A)/R;
    froerr_unstable(i) = norm(A-C*U*R,'fro');
end


semilogy(2:100,froerr_unstable,'-b','LineWidth',1);
hold on;
semilogy(2:100,froerr,'-r','LineWidth',1);
h_legend = legend('Naive CUR','StableCUR');
set(h_legend,'FontSize',20);
xlabel('number of columns and rows(c=r)','FontSize',15);
ylabel('||A-CUR||_F/||A-A_k||_F','FontSize',15);



