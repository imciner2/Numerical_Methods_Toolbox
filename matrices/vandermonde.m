function [ V ] = vandermonde( x, n, varargin )
%VANDERMONDE Compute the Vandermonde matrix of the polynomial basis
%   
% Create the Vandermonde matrix associated with the selected polynomial
% basis.
%
%
% Usage:
%   [ V ] = VANDERMONDE( x, n )
%   [ V ] = VANDERMONDE( x, n, poly)
%   [ V ] = VANDERMONDE( x, n, poly, min, max);
%
% Inputs:
%   x - Point at which to compute the polynomial
%   n - Order of the polynomial to compute
%   poly - Which polynomial basis to use for computation.
%          Possibilities include:
%              'Monomial'    - Standard monic polynomial basis (default)
%              'Chebyshev'   - 1st order Chebyshev polynomial basis
%              'Legendre'    - Legendre polynomial basis
%              'SSChebyshev' - Scaled and shifted Chebyshev polynomial
%                              basis. (Requires min and max)
%              'SSLegendre'  - Scaled and shifted Legendre polynomial
%                              basis. (Requires min and max)
%   min - The lower bound of the interval to shift the polynomial to
%   max - The upper bound of the interval to shift the polynomial to
%
% Outputs:
%   V - Vandermonde matrix
%
%
% Created by: Ian McInerney
% Created on: January 8, 2018
% Version: 1.2
% Last Modified: February 9, 2018
%
% Revision History
%   1.0 - Initial release
%   1.1 - Added SSChebyshev support
%   1.2 - Added Legendre polynomial support


%% Select the polynomial to use to build the Vandermonde matrix
poly = @monomial;
if (nargin >= 3)
    switch( varargin{1} )
        case 'SSLegendre'
            poly = @(xi, ni) sslegendrePoly(xi, ni, varargin{2}, varargin{3});
        case 'Legendre'
            poly = @legendrePoly;
        case 'SSChebyshev'
            poly = @(xi, ni) sschebyshevPoly(xi, ni, varargin{2}, varargin{3});
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

