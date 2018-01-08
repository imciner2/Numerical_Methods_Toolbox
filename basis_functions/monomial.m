function [ val, T ] = monomial( x, n )
%monomial Evaluate the monic polynomial basis at x
%   
% Compute the monomial basis of of degree n at the given point x. If x
% is multi-dimensional, then the product-monomial polynomial is computed
% and returned.
%
% Usage:
%   [ val ] = monomial( x, n )
%   [ val, T ] = monomial( x, n)
%
% Inputs:
%   x - Point at which to compute the polynomial
%   n - Order of the polynomial to compute
%
% Outputs:
%   val - Result for the highest order polynomial
%   T   - Results for all orders up-to and including the one requested
%
%
% Created by: Ian McInerney
% Created on: January 8, 2018
% Version: 1.0
% Last Modified: January 8, 2018
%
% Revision History
%   1.0 - Initial release


%% Determine the dimension of the data and create the matrix
rx = reshape(x, [], 1);     % make the points into a column vector
[d, ~] = size(rx);
T = ones(1, n+1);
inter = zeros(d, n+1);


%% Create the initial order
inter(:,1) = 1;


%% Create the remaining orders
for (i=1:1:(n))
    inter(:,i+1) = rx.^i;
end


%% Compute the column-wise product to get the final polynomial values
for (i=1:1:d)
    T = T.*inter(i,:);
end


%% The polynomial's value is the last element of the array
val = T(end);

