function [b,num_steps] = myzero_new_v2(a,b,t,f,red,options)
%[b,num_steps] = myzero_new_v2(a,b,t,@func_handle,red)
%
% Given a function handle, @func_handle, real valued scalars a and b,
% where f(a)f(b) < 0, a tolerance, t, and a reduction, red, performs 
% a modified Brent's Method. This modification guarantees O(n) performance 
% in the worst case, where n is the expected number of iterations required
% for the bisection method to converge


%Calculate intial points
fa = f(a);
fb = f(b);

if (fa<0) == (fb<0) || isnan(fa) || isnan(fb)
    error('function must have differing sign on interval points');
end

if isempty(t)
    t=eps;
end

fc = fb;
c = b;
num_steps = 1;
di = 0;
bad_count=0;
bad_max=5;

old_size=abs(b-a);

%Check options (currently only supports displaying information about
%iteratitions)
if nargin == 6
    if options.display == 1
        di = 1;
    end
end

%Display information about steps if asked to
if di == 1
    disp(' Func-count     b           f(b)           Procedure');
end
procedure = 'initial';

%Main Loop, executes while we do not have a value b where f(b) == 0 and
%the two interval points are not identical
while fb ~= 0 && a ~= b 
    %If fb and fc have the same sign then replace c with a to ensure that
    %f(b) and f(c) have opposite sign
    if (fb >0) == (fc > 0)
        c = a;      fc = fa; 
    end
    
    %If f(c) is closer to 0 than f(b) then exchange them so that f(b) is
    %always the best approximation of the root that we have found so far
    if abs(fc)<abs(fb)
        a = b;      b = c;      c = a;
        fa = fb;    fb = fc;    fc = fa;
    end
    
    %Compute tolerance and size of interval
    tol = 2.0*t*max(abs(b),1.0);
    m = 0.5*(c-b);
    
    %Check to see if interval has decreased by an appropriate amount
    if 2*abs(m)<=0.5*old_size
        old_size=2*abs(m);
        bad_count=0;
    else
        bad_count=bad_count+1;
    end
    
    %If interval size is small enough, or f(b) is identically 0, then
    %terminate
    if abs(m)<=tol || fb == 0
        break
    end
    
    %Display step information if requested
    if di == 1
        fprintf('%5.0f   %13.6g %13.6g        %s\n',num_steps,b,fb,procedure);
    end
    
    s = fb/fa;
    if a == c
        %Linear interpolation
        p = 2*m*s;      q = 1-s;
    else
        %Inverse quadratic interpolation
        q = fa/fc;
        r = fb/fc;
        p = s*(2*m*q*(q-r)-(b-a)*(r-1));
        q = (q-1)*(r-1)*(s-1);
    end
    
    if p>0, q = -q;
    else    p = -p;
    end
    
    %Check to see if interpolation step is acceptable
    if (2*p)<(3*m*q-abs(tol*q)) && (abs(fb)<red*abs(fa) || strcmp(procedure,'bisection'))...
            && bad_count<bad_max
        procedure = 'interpolation';
        d = p/q;
    else
        %Else default to bisection
        procedure = 'bisection';
        d = m;
    end

    %Choose new points
    a = b;
    fa = fb;
    
    if      abs(d) > tol,   b = b+d;
    elseif  m>0,            b = b+tol; procedure=[procedure ' bad']; %#ok<AGROW>
    else                    b = b-tol; procedure=[procedure ' bad']; %#ok<AGROW>
    end
    
    fb = f(b);
    num_steps = num_steps+1;
end
        
        
           
        
        