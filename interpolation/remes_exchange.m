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
%   [ b, E ] = REMES_EXCHANGE( n, x ):
%   [ b, E ] = REMES_EXCHANGE( n, x, poly ):
%
% Inputs:
%   n    - The order of the polynomial to use
%   x    - The datapoints to fit to (each point is a row)
%   poly - The polynomial type to use. This type must be supported by the
%          vandermonde function.
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
% Last Modified: January 248, 2018
%
% Revision History
%   1.0 - Initial release


%% If a 3rd input exists, it is the polynomial to use
poly = 'Monomial';
if (nargin == 3)
    poly = varargin{1};
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

% Extract the first points
points = x(ind, :);


%% Loop doing the algorithm until complete
STOP = 0;
while (~STOP)
    % Create the linear system matrix for the solver
    V = vandermonde(points(:,1), n+1, poly);
    alt = cos( [1:1:n+2]*pi );  % Create a vector of ones that is alternating
    V = [V, alt'];
    
    % Solve the linear system matrix
    a = V\points(:,2);
    E = a(end);
    b = a(1:end-1);
    
    % Find the residual error at each point in the set
    Vp = vandermonde(x(:,1), n+1, poly);
    P = Vp*b;
    r = P - x(:,2);
    
    % Check to see if the residual errors are all less than E
    if ( (max(abs(r)) - abs(E)) < 1e-6 )
        STOP = 1;
        break;
    end
    
    % Find the local extrema of the residual not on the ends
    for ( i = 1:1:(n+2) )
        curInd = ind(i);
        
         % Find the residual terms that are the closest ones of the same
         % sign
         ind(i) = findExtrema(r, curInd);
    end
    
    % See if there are any opposite signs outside the current end points
    % The right side
    rightPoints = (ind(end)+1):1:length(r);
    signChange = abs( sum( sign( r(rightPoints) ) ) );
    if ( signChange ~= length(rightPoints) )
        % There is another interval to the right of the end point
        % Find the largest residual in that interval
        curInd = length(r);
        tempInd = findExtrema(r, curInd);
        
        % See if that residual is larger than the one at the first point
        if ( abs(r(tempInd)) > abs(r( ind(1) )) )
            % If the new point is larger, replace the first point with it
            ind(1) = tempInd;
        end
    end
    
    % Sort the indices to make sure they are in order
    ind = sort(ind);
    
    % The left side
    leftPoints = 1:1:(ind(1));
    signChange = abs( sum( sign( r(leftPoints) ) ) );
    if ( signChange ~= length(leftPoints) )
        % There is another interval to the right of the end point
        % Find the largest residual in that interval
        curInd = 1;
        tempInd = findExtrema(r, curInd);
        
        % See if that residual is larger than the one at the first point
        if ( abs(r(tempInd)) > abs(r( ind(end) )) )
            % If the new point is larger, replace the first point with it
            ind(end) = tempInd;
        end
    end
    
    % Sort and extract the points for the next iteration
    ind = sort(ind);
    points = x(ind, :);
end


end

%% This function is designed to find the extrema of the residual in the interval
% around the current index
function [ maxInd ] = findExtrema(r, curInd)
    % Find the interval to search
    sgnInd = find( sign(r) == sign( r(curInd) ) );
    curLoc = find( sgnInd == curInd );

    % Create a list of the neighboring points
    neighList = [(curInd-curLoc+1):1:(curInd-1),... % Fill in up to the current index
                  curInd:1:(curInd + length(sgnInd(curLoc:end))-1)]';    % Fill in after the current index
    neighInd = sgnInd( sgnInd == neighList );

    % Exchange the current point in the interval with the one that has
    % the largest residual
    [~, temp] = max( abs(r(neighInd)) );
    
    % Get the index of the extrema to return
    maxInd = neighInd(temp);
end
