function [ val, T ] = legendrePoly( x, n )
%LEGENDREPOLY Evaluate the legendre polynomial at point x
%   
% Compute the Legendre polynomial of the of degree n at the
% given point x. If x is multi-dimensional, then the product-Legendre
% polynomial is computed and returned.
%
% Usage:
%   [ val ] = LEGENDREPOLY( x, n )
%   [ val, T ] = LEGENDREPOLY( x, n)
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
% see also SSLEGENDREPOLY
%
% Created by: Ian McInerney
% Created on: February 9, 2018
% Version: 1.0
% Last Modified: February 9, 2018
%
% Revision History
%   1.0 - Initial release


%% Determine the dimension of the data and create the matrix
rx = reshape(x, [], 1);     % make the points into a column vector
[d, ~] = size(rx);
T = ones(1, n+1);
inter = zeros(d, n+1);


%% Create the initial value
inter(:,1) = 1;


%% If the order is greater than 0, create the second value
if (n > 0)
    inter(:,2) = rx;
end


%% Create the remaining orders
for (i=2:1:(n))
    t = i-1;
    a = (2*t+1)/(t+1);
    b = t/(t+1);
    inter(:,i+1) = a.*rx.*inter(:,i) - b.*inter(:,i-1);
end


%% Compute the column-wise product to get the final polynomial values
for (i=1:1:d)
    T = T.*inter(i,:);
end


%% The Chebyshev polynomial's value is the last element of the array
val = T(end);

