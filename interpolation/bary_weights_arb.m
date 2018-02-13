function [ w, C ] = bary_weights_arb( x, n )
%BARY_WEIGHTS_ARB Compute the weights for barycentric interpolation using
%arbitrary interpolation points
%
% This function will compute the weights for barycentric interpolation
% using the interpolation points provided. It will automatically scale the
% weights to reduce the possibility of overflow/underflow issues in the
% floating-point computations. Information on the computation and scaling
% of the weights can be found in:
%   J.-P. Berrut and L. N. Trefethen, “Barycentric Lagrange Interpolation,”
%   SIAM Review, vol. 46, no. 3, pp. 501–517, 2004.
%
%
% Usage:
%   [ w, C ] = BARY_WEIGHTS_ARB( x, n );
%
% Inputs:
%   x - The points to interpolate through
%   n - The order of the polynomial to create (must be 1 less than the
%       number of points)
%
% Outputs:
%   w - The computed weights (in order of the points provided)
%   C - The scaling factor used to scale the difference computation
%
%
% see also BARY_COMPUTEINTERP
%
% Created by: Ian McInerney
% Created on: February 13, 2018
% Version: 1.0
% Last Modified: February 13, 2018
%
% Revision History
%   1.0 - Initial release


%% Counting starts at 0, so add 1 to n to get n+1 points
n = n+1;


%% Make sure enough points are supplied
np = length(x);
if ( np ~= n )
    error('bary_arbWeights:Mismatch between number of supplied points and polynomial degree');
end


%% Figure out the capacity of the interval
mi = min(x);
ma = max(x);
C = (ma - mi)/4;


%% Compute the weights
w = ones(n, 1);
for (j=1:1:n)
    for (k=1:1:n)
        % Skip if the points are the same points
        if (j == k)
            continue
        end
        
        di = (x(j) - x(k))/C;
        
        % Perform the multiplication
        w(j) = w(j).*di;
    end
end


%% Invert the weights to get the final result
w = 1./w;


end

