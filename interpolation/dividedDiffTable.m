function F = dividedDiffTable(x, f, n)
%DIVIDEDDIFFTABLE Calculate a divided difference table
%
% Calculate a divided difference table for a polynomial of degree n
% using points x and function values f.
%
%
% Usage:
%   [ F ] = DIVIDEDDIFFTABLE(x, f, n)
%
% Inputs:
%   x - Vector containing the x coordinates of the points
%   f - Vector containing the function values at the points
%   n - Order of the polynomial that can be fit using this table
%   
% Outputs:
%   F - Divided difference table
%
%
% Created by: Ian McInerney
% Created on: September 21, 2016
% Version: 1.0
% Last Modified: February 12, 2018
%
% Revision History
%   1.0 - Initial release


%% Add 1 to n to simplify the expressions
n = n+1;


%% Make sure the user supplied the right number of points
if (n ~= length(x))
    error('dividedDiffTable:Mismatch between number of interpolation points supplied and desired polynomial order');
end

%% Create an empty F matrix to return
F = zeros(n);


%% Fill in the first column of the matrix with function values
for ( i = 1:1:(n) )
    F(i,1) = f(i);
end


%% Fill in the rest of the matrix
for ( i = 2:1:n )
    % Iterate over the rows of F
    for ( j = 2:1:i )
        % Iterate over the columns in the row of F
        F(i,j) = (F(i, j-1) - F(i-1,j-1))/(x(i) - x(i-j+1));
    end
end

end