function [ val, T ] = sschebyshevPoly( x, n, mi, ma )
%SSCHEBYSHEVPOLY Evaluate the 1st order shifted and scaled Chebyshev
%polynomial at point x
%   
% Compute the 1D Chebyshev polynomial of the first order of degree n at the
% given point x. This Chebyshev polynomial will be shift the provided
% points such that they are contained in the range [-1,1] (where the
% Chebyshev polynomial is valid), and scale the polynomial to be 1 at the
% x=0.
%
%
% Usage:
%   [ val, T ] = SSCHEBYSHEVPOLY( x, n, mi, ma )
%
% Inputs:
%   x  - Point at which to compute the polynomial
%   n  - Order of the polynomial to compute
%   mi - The minimum of the interval
%   ma - The maximum of the interval
%
% Outputs:
%   val - Result for the highest order polynomial
%   T   - Results for all orders up-to and including the one requested
%
%
% Created by: Ian McInerney
% Created on: January 29, 2018
% Version: 1.0
% Last Modified: January 29, 2018
%
% Revision History
%   1.0 - Initial release


%% Determine the dimension of the data and create the matrix
rx = reshape(x, [], 1);     % make the points into a column vector
d = length(x);
inter_num = zeros(d, n+1);
inter_den = zeros(d, n+1);


%% Create the initial value
inter_num(:,1) = ones(d,1);
inter_den(:,1) = ones(d,1);


%% Create the constants for shifting and scaling
al = (ma + mi)/2;
rh = (ma - mi)/2;


%% Shift and scale x to create the numerator terms
sx = (rx - al)./rh;
de = -al/rh;


%% If the order is greater than 0, create the second value
if (n > 0)
    inter_num(:,2) = sx;
    inter_den(:,2) = de;
end


%% Create the remaining orders
for (i=2:1:(n))
    inter_num(:,i+1) = 2.*sx.*inter_num(:,i) - inter_num(:,i-1);
    inter_den(:,i+1) = 2.*de.*inter_den(:,i) - inter_den(:,i-1);
end


%% Compute the final polynomial values
T = inter_num./inter_den;


%% The Chebyshev polynomial's value is the last element of the array
val = T(:,end);

