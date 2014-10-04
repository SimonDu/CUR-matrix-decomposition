function Afft = realfft(A)
% Afft = realfft(A)
%
% Returns the normalized real fft applied to the columns of A, a matrix 
% with n rows
% Assume n>=5 say, so we don't have to worry about edge cases

% Based on the observation
% If f is real, then f = fft(x) satisfies
% f(1) = sum(x)
% let n = length(x)
% if n odd, then let a = 2:ceil(n/2) and b = ceil(n/2)+1 : n
%    then f(a) = conjugate(rev(b))
% if n even, then let a = 2:n/2 and b=n/2+2:n
%    then f(a) = conjugate(rev(b)) and f(n/2+1) = sum(x(i)*(-1)^{i-1},
%    i=1,..n)

    n = size(A,1);
    Afft = 1/sqrt(n)*fft(A);
    if rem(n,2) == 1
        realindices = 2:ceil(n/2);
        imagindices = (ceil(n/2)+1):n;
        Afft(realindices, :) = sqrt(2)*real(Afft(realindices, :));
        Afft(imagindices, :) = sqrt(2)*imag(Afft(imagindices, :));
    else
        realindices = 2:(n/2);
        imagindices = (n/2+2):n;
        Afft(realindices, :) = sqrt(2)*real(Afft(realindices, :));
        Afft(imagindices, :) = sqrt(2)*imag(Afft(imagindices, :));
    end
    
end
