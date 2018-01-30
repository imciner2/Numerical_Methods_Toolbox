function [ b, E ] = remes_exchange( n, x, varargin )
%REMES_EXCHANGE Remes Exchange algorithm to compute the minimax polynomial
%over a discrete set
%
% Find the minimax polynomial using a discrete set of points as the
% function to match instead of a continuous function. This function uses
% the Remes Exchange Algorithm to compute the polynomial. An overview of
% this method and the minimax problem can be found inside:
%   W. Fraser, “A Survey of Methods of Computing Minimax and Near-Minimax
%   Polynomial Approximations for Functions of a Single Independent
%   Variable,” J. ACM, vol. 12, no. 3, pp. 295–314, 1965.
%
%
% Usage:
%   [ b, E ] = REMES_EXCHANGE( n, x );
%   [ b, E ] = REMES_EXCHANGE( n, x, poly );
%   [ b, E ] = REMES_EXCHANGE( n, x, poly, MITER );
%   [ b, E ] = REMES_EXCHANGE( n, x, poly, MITER, origin );
%
% Inputs:
%   n      - The order of the polynomial to use
%   x      - The datapoints to fit to (each point is a row)
%   poly   - The polynomial type to use. This type must be supported by the
%            vandermonde function (default Monomial).
%   MITER  - The maximum number of iterations to use (default 1000).
%   origin - Force this value at the origin
%
% Outputs:
%   b - The polynomial coefficients (lowest ordered first)
%   E - The minimax error for the polynomial
%
% See also VANDERMONDE.
%
% Created by: Ian McInerney
% Created on: January 24, 2018
% Version: 1.0
% Last Modified: January 29, 2018
%
% Revision History
%   1.0 - Initial release


%% Parse the input
p = inputParser;
addOptional(p, 'poly', 'Monomial', @ischar);
addOptional(p, 'MAX_ITER', 1000);
addOptional(p, 'origin', NaN);
parse(p, varargin{:});

poly = p.Results.poly;
MAX_ITER = p.Results.MAX_ITER;
origin = p.Results.origin;

if ( strcmp(poly, 'SSChebyshev') == 1)
    van = @(x) vandermonde(x, n+1, poly, min(x(:,1)), max(x(:,1)));
else
    van = @(x) vandermonde(x, n+1, poly);
end



%% Determine how many points were provided and make sure there were enough
[m, ~] = size(x);
if ( (n+2) > m )
    error('remes_exchange: Not enough points provided for the desired polynomial order');
end


%% Make sure that the points are sorted in ascending order
x = sortrows(x, 1);


%% Create the first set of points to use for the approximation
% Take as equally spaced points as possible in the interval
ind = linspace(1, m, n+2);
ind = round(ind);
oldInd = ind;

% Extract the first points
points = x(ind, :);
np = length(ind);


%% Loop doing the algorithm until complete
STOP = 0;
k = 1;
while (~STOP)
    % Create the linear system matrix for the solver
    V = van( points(:,1) );
    alt = cos( [1:1:np]*pi );  % Create a vector of ones that is alternating
    V = [V, alt'];
    
    % Setup the linear system
    if ( isnan(origin) )
        % Let the constant parameter be free
        c = points(:,2);
        b = [];
    else
        % Force a value at the origin of the polynomial
        % This is equivalent to fixing the constant parameter and removing
        % it from the other side
        V = V(:, 2:end);
        c = points(:,2) - origin.*ones(np, 1);
        b = origin;
    end
    
    % Solve the linear system and create the coefficient vector
    a = V\c;
    b = [b; a(1:end-1)];
    
    % This is the approximation error for this try
    E = abs( a(end) );
    
    % Find the residual error at each point in the set
    Vp = van( x(:,1) );
    P = Vp*b;
    r = P - x(:,2);
    
    % Check to see if the residual errors are all less than E
    if ( (max(abs(r)) - E) < 1e-4 )
        STOP = 1;
        break;
    end
    
    ind = findExtrema(r);
        
    % Sort and extract the points for the next iteration
    ind = sort(ind);
    points = x(ind, :);
    np = length(ind);
    
    % See if the maximium iteratons has been reached
    if (k == MAX_ITER)
        warning('remes_exchange:Reached maximum iteration limit');
        STOP = 1;
        break;
    end
    k = k + 1;
    
    % See if the indices have not changed between last time and this time
    if ( length(oldInd) == length(ind) )
        if (oldInd == ind)
            warning('remes_exchange:Reached stationary extrema set without meeting tolerance');
            STOP = 1;
            break;
        end
    end
    
    % Save the new indices for comparison
    oldInd = ind;
end


end

%% This function is designed to find the extrema of the residual in the interval
% around the current index
function [ ind ] = findExtrema(r)
    % Find the intervals to search
    sig = sign(r);
    intEnd = find( sig(1:end-1) ~= sig(2:end) );
    intStart = intEnd + 1;
    
    % Add the start and end of the array to the values
    intEnd = [intEnd; length(r)];
    intStart = [1; intStart];
    
    % Iterate over every interval that was found
    for ( i=1:1:length(intStart) )
        % Pull out the residuals
        rInt = r( intStart(i):1:intEnd(i) );
        
        % Find the largest one and get its index
        [~, mr] = max( abs(rInt) );
        ind(i) = mr + intStart(i) - 1;
    end
end
