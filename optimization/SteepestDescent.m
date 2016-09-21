function [ x, NumIters ] = SteepestDescent( fHandle, dFhandle, x0, TOL, MaxIters )
%STEEPESTDESCENT Perform a steepest descent optimization on the given
%function
%
%   This function will perform a steepest descent optimization on the given
%   unconstrained function
%
%   This function takes 5 arguments:
%       fHandle - The function handle to the function to compute f(x)
%       dFhandle - The function hanle to the function to compute the
%                  gradiant of f(x)
%       x0 - The inital guess
%       TOL - The tolerance to solve to
%       MaxIters - The maximum number of points to search
%
%   This function returns 2 items:
%       x - The optimum point found
%       NumIters - The number of iterations used
%
%   Author: Ian McInerney

%% Check to make sure inputs make sense
if (MaxIters < 0)
    error('SteepestDescent::Max Iterations must be positive');
end


%% Setup initial variables
x = x0;
mStop = 0;
NumIters = 1;


%% Run the algorithm until a specified tolerance is reached
while ( (NumIters < MaxIters) && ~mStop )
    % Find the direction to search in
    p = -1*feval(dFhandle, x);
    
    % Perform a backtracking line search to get the step length
    alpha = LineSearch_Backtracking(fHandle, dFhandle, p, x, 1000);
    
    % Compute the next point
    xOld = x;
    x = xOld + alpha*p;
    
    % Check the tolerance to see if we can stop the algorithm
    if ( norm(x - xOld, 2) < TOL)
        mStop = 1;
    end
end


%% Check if the algorithm timed out, warn if it did
if (NumIters == MaxIters)
    warning('SteepestDescent::Maximum Number of Iterations Reached');
end

end

