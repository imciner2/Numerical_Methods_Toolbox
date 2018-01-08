function [ V ] = vandermonde( x, n, varargin )
%VANDERMONDE Compute the Vandermonde matrix of the polynomial basis
%   
% Create the Vandermonde matrix associated with the selected polynomial
% basis.
%
%
% Usage:
%   [ V ] = vandermonde( x, n )
%   [ V ] = vandermonde( x, n, poly)
%
% Inputs:
%   x - Point at which to compute the polynomial
%   n - Order of the polynomial to compute
%   poly - Which polynomial basis to use for computation.
%          Possibilities include:
%              'Monomial' - Standard monic polynomial basis (default)
%              'Chebyshev' - 1st order Chebyshev polynomial basis
%
% Outputs:
%   V - Vandermonde matrix
%
%
% Created by: Ian McInerney
% Created on: January 8, 2018
% Version: 1.0
% Last Modified: January 8, 2018
%
% Revision History
%   1.0 - Initial release


%% Select the polynomial to use to build the Vandermonde matrix
poly = @monomial;
if (nargin == 3)
    switch( varargin{1} )
        case 'Chebyshev'
            poly = @chebyshevPoly;
        case 'Monomial'
            poly = @monomial;
        otherwise
            poly = @monomial;
    end        
end


%% Create the matrix structure
[numPoints, ~] = size(x);
V = zeros(numPoints, n);


%% Construct the Vandermonde matrix
for ( i=1:1:numPoints )
    [ ~, V(i,:) ]= poly(x(i,:), n-1);
end

end

