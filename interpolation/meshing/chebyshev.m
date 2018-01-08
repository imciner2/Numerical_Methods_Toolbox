function [ chebyNodes ] = chebyshev( n, lb, ub, varargin )
%CHEBYSHEV Create the Chebyshev nodes
%
% Computes the Chebyshev nodes associated with the Chebyshev polynomial
% of the first order. For multiple dimensions, the roots are computed in
% a tensorized way.
%
%
% Usage:
%   [ chebyNodes ] = chebyshev( n, lb, ub )
%   [ chebyNodes ] = chebyshev( n, lb, ub, d )
%
% Inputs:
%   n  - The number of roots to compute
%   lb - Lower bound for the roots
%   ub - Upper bound for the roots
%   d  - The number of dimensions to create roots in
%
%   For multiple dimensions, lb and ub can be a vector with the respective
%   bounds for the ith dimension specified by the ith element.
%
% Outputs:
%   chebyNodes - Array containing the computed nodes
%
% Created by: Ian McInerney
% Created on: January 7, 2018
% Version: 1.0
% Last Modified: January 7, 2018
%
% Revision History:
%   1.0 - Initial release

%% the number of dimensions to create the Chebyshev nodes in
d = 1;
if (nargin == 4)
    d = varargin{1};
end


%% Create the bounds array (if needed)
if ( length(lb) ~= d )
    lb = lb(1).*ones(1,d);
end
if ( length(ub) ~= d )
    ub = ub(1).*ones(1,d);
end


%% Create the Chebyshev nodes

% Loop to create the raw points
origCheby = zeros(n,1);
for (k=1:1:n)
    origCheby(k) = cos((2*k - 1)/(2*n)*pi);
end

%% Loop over each dimension moving the points
chebyNodes = [];
for( i = 1:1:d )
    % Create the Chebyshev points for the dimension
    bias(1) = (ub(i) - lb(i))/2;
    bias(2) = (ub(i) + lb(i))/2;
    cheby = bias(2).*ones(n,1) + bias(1).*origCheby;
    
    % Iterate over each point and replicate it accordingly
    temp = [];
    for (k=1:1:n)
        temp = [temp; repelem(cheby(k), n^(i-1))'];
    end
    
    % Replicate the previous points to fill the dimension
    chebyNodes = repmat(chebyNodes, n, 1);
    
    % Concatenate the repated points with the new points
    chebyNodes = [chebyNodes, temp];
end


end

