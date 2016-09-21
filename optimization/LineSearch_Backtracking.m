function [ alpha ] = LineSearch_Backtracking( fHandle, dFhandle, p, x, MaxIters )
%LINESEARCH_BACKTRACKING Perform a backtracking line search
%   This function performs a backtracking line search to find a value of
%   alpha that decreases the value of F in the search direction such that
%   the sufficient decrease condition is satisfied.
%
%   This function takes 5 arguments:
%       fHandle - The function handle to the function to compute f(x)
%       dFhandle - The function hanle to the function to compute the
%                  gradiant of f(x)
%       p - The unit-vector giving the search direction
%       x - The starting point for the search
%       MaxIters - The maximum number of points to search
%
%   This function returns 1 argument:
%       alpha - The alpha found to have satisfied the sufficient decrease
%               condition
%
%   Author: Ian McInerney


%% Setup some variables necessary for the algorithm

% These are some simple constants for the algorithm
c = 1E-4;       % Sufficient decrease constant
rho = 0.3;      % Step variation constant

% Find the initial function value and initial step length
Fstart = feval(fHandle, x);
gradF = feval(dFhandle, x);
alpha = 0.1;


% Create the first sufficient condition
sufficientCondition = Fstart + c*alpha*gradF'*p;
i = 1;


%% Run the backtracking algorithm to find the step length
% This loop will run until the function value is less than the desired
% function value (meets the sufficient decrease condition)
F = Fstart;
while ( F >= sufficientCondition && i < MaxIters)
    % Update alpha using the rho constant
    alpha = rho*alpha;
    
    % Recompute the function with the new step length and distance
    F = feval(fHandle, x + alpha.*p);
    
    % Recompute the sufficient condition with the new alpha
    sufficientCondition = Fstart + c*alpha*gradF'*p;
    
    % Increment the iteration counter
    i = i + 1;
end

if (i == MaxIters)
    warning('LineSearch_Backtracking::Maximum Iteration Limit Reached');
end

end
