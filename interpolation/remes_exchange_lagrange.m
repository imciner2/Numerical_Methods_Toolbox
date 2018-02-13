function [ b, E ] = remes_exchange_lagrange( n, x, varargin )
%REMES_EXCHANGE_LAGRANGE Remes Exchange algorithm to compute the minimax
% interpolation polynomial of the Lagrange type over a discrete set
%
% Find the minimax polynomial using a discrete set of points as the
% function to match instead of a continuous function. This function uses
% the Remes Exchange Algorithm to compute the polynomial. An overview of
% this method and the minimax problem can be found inside:
%   W. Fraser, “A Survey of Methods of Computing Minimax and Near-Minimax
%   Polynomial Approximations for Functions of a Single Independent
%   Variable,” J. ACM, vol. 12, no. 3, pp. 295–314, 1965.
%
% This function uses the barycentric method of evaluating the interpolating
% polynomial over the arbitrary set of points provided. This method is used
% to reduce numerical instability in the polynomial computation.
%
%
% Usage:
%   [ b, E ] = REMES_EXCHANGE_LAGRANGE( n, x );
%   [ b, E ] = REMES_EXCHANGE_LAGRANGE( n, x, MITER );
%   [ b, E ] = REMES_EXCHANGE_LAGRANGE( n, x, MITER, RITER );
%   [ b, E ] = REMES_EXCHANGE_LAGRANGE( n, x, MITER, RITER, TOL );
%
% Inputs:
%   n     - The order of the polynomial to use
%   x     - The datapoints to fit to (each point is a row)
%   MITER - The maximum number of iterations to use (default 1000)
%   RITER - The number of iterations to run before restarting (default
%           500). To disable restarting, pass in 0.
%   TOl   - The relative tolerance for ending the algorithm. The algorithm
%           will terminate once the change in error is below this value
%           (default 1e-5).
%
% Outputs:
%   b - The points used in the minimax polynomial
%   E - The minimax error for the polynomial
%
% See also REMES_EXCHANGE.
%
% Created by: Ian McInerney
% Created on: February 12, 2018
% Version: 1.0
% Last Modified: February 12, 2018
%
% Revision History
%   1.0 - Initial release


%% Parse the input
p = inputParser;
addOptional(p, 'MAX_ITER', 1000);
addOptional(p, 'RES_ITER', 500);
addOptional(p, 'TOL', 1e-5);
parse(p, varargin{:});

MAX_ITER = p.Results.MAX_ITER;
RES_ITER = p.Results.RES_ITER;
TOL      = p.Results.TOL;


%% Determine how many points were provided and make sure there were enough
[m, ~] = size(x);
if ( (n+1) > m )
    error('remes_exchange_lagrange: Not enough points provided for the desired polynomial order');
end


%% Make sure that the points are sorted in ascending order
x = sortrows(x, 1);

E = NaN;

%% Create the first set of points to use for the approximation
% Take as equally spaced points as possible in the interval
ind = linspace(1, m, n+1);
ind = round(ind);
oldInd = ind;

% Extract the first points
points = x(ind, 1);
np = length(ind);


%% Loop doing the algorithm until complete
STOP = 0;
k = 1;
while (~STOP)
    % Create the barycentric weights for the current point set
    [w, C] = bary_weights_arb(points, np-1);

    % The error function has alternating sign of magnitude 1 at each extrema
    f = cos(pi*[0:1:np-1]);
    
    % Compute the polynomial value at each x location and scale it
    interp_den = bary_computeInterp(0, points, f, w, C);
    interp_num = bary_computeInterp(x(:,1), points, f, w, C)./interp_den;
    
    % Find the residual error at each point in the set
    r = interp_num - x(:,2);
    newE = max( abs(r) );
    
    % Stop once E is no longer decreasing very much
    if ( abs(newE - E) < TOL )
        STOP = 1;
        break;
    end
    E = newE;
    
    % Find the new extrema in the residual set to update the indices
    newInd = findExtrema(r);
 
    % See if too many indices are being returned
    numInd = length(newInd);
    while ( numInd > np )
        % Choose a random index to remove
        remInd = round( (numInd-1)*rand(1)+1 );
        
        % Create a logical indexer that ignores that index
        remArr = ones(numInd, 1);
        remArr(remInd) = 0;
        
        % Remove the index
        newInd = newInd( logical(remArr) );
        numInd = length(newInd);
    end
    
    % See if too few indices are being returned
    while ( length(newInd) < np )
        % If there are not, then choose another point at random to use
        temp = round( (m-1)*rand(1)+1 );
        newInd = unique( [newInd temp] );
    end
    
    % Sort and extract the points for the next iteration
    ind = sort(newInd);
    points = x(ind, 1);
    np = length(ind);
    
    % See if the maximium iteratons has been reached
    if (k == MAX_ITER)
        warning('remes_exchange_lagrange:Reached maximum iteration limit');
        STOP = 1;
        break;
    end
    k = k + 1;
    
    % Randomly choose new points if restarting is requested
    if ( mod(k, RES_ITER) == 0 )
        ind = randperm(m, np);
        ind = sort(ind);
        points = x(ind, 1);
    end
    
    % See if the indices have not changed between last time and this time
    if ( length(oldInd) == length(ind) )
        if (oldInd == ind)
            STOP = 1;
            break;
        end
    end
    
    % Save the new indices for comparison
    oldInd = ind;
end

% Save the points used for interpolating to pass back to the caller
b = points;
E = max( abs(r) );

end

%% This function is designed to find the extrema of the residual
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
