function [ p ] = bary_computeInterp( xbar, x, f, w, varargin )
%BARY_COMPUTEINTERP Compute the polynomial interpolation using the
%Barycentric formula
%
% This function will perform function interpolation using Lagrange
% polynomials through the barycentric interpolation formula of the second
% form.
%
% Interpolation using this method is usually more numerically stable than
% the Divided Difference/Horner's method for the Newton polynomial
% evaluation.
%
% The weights used in this interpolation are pre-computed using one of the
% other functions.
%
% For more about the barycentric interpolation method, see:
%   [1] J.-P. Berrut and L. N. Trefethen, “Barycentric Lagrange Interpolation,”
%       SIAM Review, vol. 46, no. 3, pp. 501–517, 2004.
%   [2] N. J. Higham, “The numerical stability of barycentric Lagrange
%       interpolation,” IMA Journal of Numerical Analysis, vol. 24, no. 4,
%       pp. 547–556, 2004.
%
%
% Usage:
%   [ p ] = BARY_COMPUTEINTERP( xbar, x, f, w );
%   [ p ] = BARY_COMPUTEINTERP( xbar, x, f, w, C );
%
% Inputs:
%   xbar - The points to compute the interpolant at
%   x    - The points where the interpolation function passes through
%   f    - The function value at the x points
%   w    - The weights
%   C    - Scaling of the difference (necessary if weight computation used
%          it).
%
% Outputs:
%   p - The interpolated polynomial values
%
%
% see also BARY_WEIGHTS_ARB
%
% Created by: Ian McInerney
% Created on: February 13, 2018
% Version: 1.0
% Last Modified: February 13, 2018
%
% Revision History
%   1.0 - Initial release

%% Make sure enough interpolation points are supplied
n = length(w);
if ( length(x) ~= n )
    error('bary_computeInterp:Mismatch between interpolation points and weights');
end


%% Parse the inputs
p = inputParser;
addOptional(p, 'C', 1);
parse(p, varargin{:});

C = p.Results.C;


%% Iterate over the points to get the interpolation at each one
np = length(xbar);

% Initialize the variables
num = zeros(np, 1);
den = zeros(np, 1);
exa = zeros(np,1);

% Iterate over the weights
for (j=1:1:n)
    xdi = (xbar - x(j))./C;
    coe = w(j)./xdi;
    num = num + coe.*f(j);
    den = den + coe;
    exa( xdi==0 ) = j;
end


%% Divide to get the final result
p = num./den;


%% Fix places where NaN occur (e.g. if x=x(j))
j = find(exa);
p(j) = f( exa(j) );

end
