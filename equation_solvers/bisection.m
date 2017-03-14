function [xSol, NumIters, x] = bisection(f, a, b, TOL, MaxIters)
%BISECTION Find a root of f using the bisection method
%
% Find a root of the function f using the bisection method
%
% Input Arguments:
%   f - Function pointer to function f(x) to compute the value of f
%   a - Left-side of search interval
%   b - Right-side of search interval
%   TOL - Stopping tolerance
%   MaxIters - The maximum number of iterations to perform
%
% Output Variables:
%   x - Sequence of guesses for the root
%   NumIters - Number of iterations before convergence occured

% Perform some input validation
fa = feval(f, a);
fb = feval(f, b);

if (a >= b)
    error('Bisection: a must be less than b');
end
if (fa*fb > 0)
    error('Bisection: No Root In Interval');
end

% Preallocate the array
x = nan(MaxIters,1);

% Start the iteration counter at the first iteration
NumIters = 0;

% Stopping variable
mstop = 0;

% Iterate and find the root
while( (NumIters < MaxIters) && ~mstop )
    NumIters = NumIters + 1;
    
    % Create the next root guess
    x(NumIters) = 0.5*(a+b);
    
    fa = feval(f, a);
    fx = feval(f, x(NumIters));
    % Check the stopping criteria, and set the flag if met
    if ( ( 0.5*(b-a) < TOL) && ( fx < TOL) )
        mstop = 1;
    elseif ( fx*fa < 0 )
        % Root is in left-half of interval, replace b
        b = x(NumIters);
    else
        % Root is in right-half of interval, replace a
        a = x(NumIters);
    end
end

if ( NumIters == MaxIters)
    warning('Bisection: Maximum Number of Iterations Reached');
end

% Truncate the array being returned to remove unused parts
x = x(1:NumIters);
xSol = x(end);

end
