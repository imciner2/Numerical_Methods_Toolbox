function [ sigma ] = svd_est_johnson( A, varargin )
%SVD_SMALLEST_JOHNSON Provide a lower bound for the Singular Values
%
% Creates a lower bound for the singular values using the estimate
% provided in:
%   Johnson, C. R. (1989). A Gersgorin-type lower bound for the smallest
%   singular value. Linear Algebra and Its Applications, 112, 1–7.
%   https://doi.org/10.1016/0024-3795(89)90583-1
%
% Usage:
%   [ sigma ] = svd_est_johnson( A )
%   [ sigma ] = svd_est_johnson( A, sat )
%
% Inputs:
%   A - The matrix to estimate the lower bound of
%   sat - Flag saying whether or not to constrain the singular value to be
%         positive (e.g. if sigma < 0, set sigma = 0)
%
% Outputs:
%   sigma - Lower bound for the smallest singular value
%
%
% Created by: Ian McInerney
% Created on: January 8, 2018
% Version: 1.0
% Last Modified: January 8, 2018
%
% Revision History:
%   1.0 - Initial release

% Get the size of A and find the smallest dimension (then make it the rows)
[n, m] = size(A);
if ( m < n )
    A = A';
    [n, m] = size(A);
end

saturate = 1;
if (nargin == 2)
    saturate = varargin{1};
end

A = A(1:n, 1:n);

% Extract the diagonal of A
vecdiagA = diag(A);
nodiagA = A - diag(vecdiagA);

% Sum along the rows and the columns of the matrix
r = sum( abs(nodiagA), 1);
c = sum( abs(nodiagA), 2);

% Average the rows and column sums
av = 0.5.*(r' + c);

% Subtract the average from the diagonal elements
res = abs( vecdiagA ) - av;

% Get the smallest value, and saturate to 0 if desired
sigma = min(res);
if (sigma < 0 && saturate == 1)
    sigma = 0;
end

end

