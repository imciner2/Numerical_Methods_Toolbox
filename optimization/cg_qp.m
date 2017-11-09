function [ x_optimal, res ] = cg_qp(A, b, x0, varargin)
%cg_qp Minimize the given quadratic function using the conjugate gradient method
%
% Minimize the unconstrained quadratic program 
%         min 1/2 x'Ax - b'x
% using the Conjugate Gradient method.
%
% This function is based on the algorithm given by Nocedal & Wright in
% their book 'Numerical Optimization', second edition.
%
% Outputs:
%   x_final - The final (optimal) x iterate
%   res - The L2-norm of the residual at each iteration
%
% Inputs:
%   A - The Hessian of the quadratic program
%   b - The coefficients for the linear term
%   x0 - The initial point
%
% Options:
%   'MaxIterations' - Specify the maximum number of iterations to perform
%                     Since this uses the Conjugate Gradient method, the
%                     number of iterations must be less than or equal to
%                     the number of variables in the problem.
%                     The default is the number of variables.
%
%   'ToleranceFunction' - Specify a custom function to use to compute if
%                         the solution is within tolerance. It must have
%                         arguments (x, r) and return a boolean (true if
%                         within tolerance, false if not).
%                         The default tolerance function uses the residual:
%                           @(x, r) ( norm( r, 2 ) < TOL )
%
%   'Tolerance' - Specify the tolerance values for the default tolerance
%                 function (if it is used).
%                 The default is to use 1e-4.
%
%
% Created by: Ian McInerney
% Created on: November 3, 2017
% Version: 1.0
% Last Modified: November 3, 2017
%
% Revision History:
%   1.0 - Initial Release

% Make sure the number of arguments is correct
numArgs = length(varargin);
if ( mod(numArgs, 2) ~= 0 )
    error('Incorrect number of arguments');
end

% Determine the number of variables in the problem
numVars = length(x0);

% Default tolerance and number of iterations
MAX_ITER = numVars;
TOL = 1e-4;

% Parse the arguments
if (numArgs > 0)
    for i=1:2:numArgs
        switch( varargin{i} )
        case 'MaxIterations'
            MAX_ITER = varargin{i+1};
            if (MAX_ITER > (numVars+1))
                warning('Requested more iterations than variables, ignoring request');
                MAX_ITER = numVars;
            end
            
        case 'Tolerance'
            TOL = varargin{i+1};

        case 'ToleranceFunction'
            tolFunc = varargin{i+1};
        end
    end 
end

% Create the tolerance function if it doesn't exist
if ( ~exist('tolFunc', 'var') )
    tolFunc = @(x, r) ( norm( r, 2 ) < TOL );
end

% Pre-allocate the arrays
res = zeros(1, MAX_ITER);

% Setup the initial variables
x = x0;
r = A*x - b;
p = -r;

res(1) = norm( r, 2);

k = 1;

%% Do the algorithm loop
STOP = 0;
while ~STOP
    % Compute the step length
    alpha = -( r'*p ) / ( p'*A*p );
 
    % Compute the next point and residual
    x_n = x + alpha*p;
    r_n = A*x_n - b;
    
    % Compute the next direction
    beta = ( r_n'*A*p ) / ( p'*A*p );
    p_n = -r_n + beta*p;
   
    % Save the residual distance
    res(k+1) = norm( r, 2 );
    
    % Save the new variables for the next loop
    x = x_n;
    r = r_n;
    p = p_n;

    % Update the iteration counter
    k = k + 1;
    
    % Stop after the max number of iterations
    % k starts at 1, then should go for MAX_ITER, stopping once it is hit
    if ( k == (MAX_ITER+2) )
        warning('Maximum interations reached before tolerance');
        STOP = 1;
    end
    
    % Check the tolerance, and if met then stop
    if ( tolFunc(x, r ) )
        STOP = 1;
    end
end


%% Output the final point
res = res(1:k);
x_optimal = x;

end