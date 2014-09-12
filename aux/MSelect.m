function pi = MSelect(U,k,l)
%Rank k approximation with Ming's algorithm
%   deterministic, unweighted algorithm for use in sparse PCA CUR
%
%   inputs:
%       U       orthogonal matrix to select cols
%               typically is V_k of rank k TSVD
%       k       rank of approximation (k<<n)
%       l       max rank of C and R (k<l<<n)
%
%   outputs:      
%       pi      index set |pi|=l s.t.
%               sigma_min(U(:,id))>=sqrt(l)-sqrt(k))/sqrt(k)

n=size(U,1);

% fixed inverse trace constant
% OLD CONSTANT:
T_old=(sqrt(k)/(sqrt(l)-sqrt(k)))*n;
% NEW CONSTANT:
T_hat=k*(n+(l+1)/2-k)+sqrt(k*l*(n-(l-1)/2)*(n+(l+1)/2-k));
T_hat=T_hat/(l-k);
F=@(t) (1-k/t)*l/(n-(l-1)/2-k+t)-k/t;
T=T_hat*(1+F(T_hat));
clear T_hat F


% rank k T-SVD
%[U0,S,U]=svds(X,k);

% initialize vars
A=sparse(k,k);
I=[];
lambda=-Inf;
lambda_hat=-Inf;
% main loop to build C and R
for r=0:l-1
    % step 1: compute lambda so tr(A-lambda*I)^-1=T
    s=eig(A);
    
    nn = 40;
    tol = 2*eps;
    fun=@(x) sum(1./(s-x))-T;
    a=max(lambda_hat,-10); % -10
    b=min(s)-(1e-10);
    if fun(a)>0 & fun(b)>0
        [lambda_new,num_steps] = myzero_new_v2(lambda,b,tol,fun,.5);
        if lambda_new-lambda_hat<-1e-6
            error('decreasing bound');
        else
            lambda=lambda_new;
        end
    end
    [lambda,num_steps] = myzero_new_v2(a,b,tol,fun,.5);
    
    % old method with 'roots'
    %{    
    TR=zeros(k+1,k+1);
    TR(1:end-1,end)=1;
    TR(end,end)=-T;
    for i=1:k
        for j=1:k-1
            ind=mod(i+j-2,k)+1;
            TR(i,:)=conv(TR(i,:),[-1,s(ind)],'same');
        end
        TR(end,:)=conv(TR(end,:),[-1,s(i)],'same');
    end
    P=sum(TR,1);
    R=real(roots(P));
    %R(R<=lambda)=inf; %%%% CHECK SELECTION CRITERIA
    lambda=min(R);
    %}


    % step 2: solve for lambda_hat
    c=n-r+sum((1-s)./(s-lambda));
    lh_update=@(lh) (lh-lambda)*c- ...
        sum((1-s)./((s-lambda).*(s-lh)))/ ...
        sum((1)./((s-lambda).*(s-lh)));
    l0=.5*(lambda+s(k));
    % OLD SOLVE
    %lambda_hat=fzero(lh_update,l0);
    % NEW SOLVE
    [lambda_hat,num_steps] = myzero_new_v2(lambda,min(s)-(1e-8),tol,lh_update,.5);
    
    
    %{
% use fzero here instead
% initial guess (lambda+lambda_k)/2 -> value should be in between
% fzero gives 1 value.  roots gives many values and does a bad job
    TR=zeros(k,k+1);
    RR=zeros(k,k+1);
    c=n-r+sum((1-s)./(s-lambda));
    %TR(1:end-1,end-1:end)=[c*ones(k-1,1),-lambda*c*ones(k-1,1)];
    %TR(end,end)=1;
    TR(:,end-1:end)=[c*ones(k,1),-lambda*c*ones(k,1)];
    for i=1:k
        RR(i,end)=1-s(i);
        for j=1:k-1
            ind=mod(i+j-1,k)+1;
            cj=s(ind)-lambda;
            TR(i,:)=conv(TR(i,:),[-cj,cj*s(ind)],'same');
            RR(i,:)=conv(RR(i,:),[-cj,cj*s(ind)],'same');
        end
    end
    P=sum(TR,1)-sum(RR,1);
    R=roots(P);
    R(R<=lambda)=-inf;%inf;;
    lambda_hat=max(R);%min(R);  %CHECK IF MAX OK
                                %an eigenvalue could be the root otherwise
%}
    %{
    % check:
    (lambda_hat-lambda)*(n-r+sum((1-s)./(s-lambda)))                           
    an=sum((1-s)./((s-lambda).*(s-lambda_hat)))                            
    ad=sum((1)./((s-lambda).*(s-lambda_hat)))
    an/ad
    %}

    % step 3: select i and update
    trinv=sum(1./(s-lambda));
    trhinv=sum(1./(s-lambda_hat));
    IC=setdiff([1:n],I);
    for ix=IC
        ui=U(ix,:)'; % transpose of U/V
        d=1+ui'*((A-lambda_hat*eye(k))\ui);
        adj=ui'*((A-lambda_hat*eye(k))\((A-lambda_hat*eye(k))\ui))/d;
        if trhinv-adj<=trinv
            A=A+ui*ui';
            I=[I,ix];
            break;
        end
    end
end

%C=U0*S*U(I,:)';
pi=I;

end

