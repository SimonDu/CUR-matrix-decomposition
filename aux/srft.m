function C = srft(A,m,c)

% my version of srft implementation...
%
% d = complex(randn(1,n), randn(1,n));
% d = bsxfun(@rdivide, d, sqrt(sum(d.*conj(d), 2)));
% D = diag(d);
% 
% F = fft(l/sqrt(n)*eye(n))/l;
% norm(F)
% 
% S = datasample(eye(n),l,2,'Replace',true);
% 
% out = D*(F*S);


%gitten's version of srft
% input
% -m, number of rows in A
colindices = randperm(m);
colindices = colindices(1:c);
d = 1-2*(rand(1,m) > .5);

% compute A*D*F' "efficiently" using realfft, the fact that A is symmetric, 
% and that A*D is 10 or so times slower than using bsxfun
% my implementation of the realfft is about 5 times slower than the regular
% fft

Y = realfft(full(bsxfun(@times,A,d)'));
C = Y(colindices, :)';