function [ coeff, E ] = lp_minimaxPoly( n, x, varargin )
%LP_MINIMAXPOLY Compute the minimax polynomial using a LP
%
% Compute the minimax polynomial approximation of order n to a discrete set
% of points. This method uses the linear programming formulation of the
% problem, and solves it using MATLAB's linprog function (requires the
% Optimization Toolbox).
%
%
% Usage:
%   [ coeff, E ] = lp_minimaxPoly( n, x );
%   [ coeff, E ] = lp_minimaxPoly( n, x, poly );
%   [ coeff, E ] = lp_minimaxPoly( n, x, poly, origin );
%   [ coeff, E ] = lp_minimaxPoly( n, x, poly, origin, alg );
%
% Inputs:
%   n      - The order of the polynomial to fit with
%   x      - The points to use to do the fitting
%   poly   - The polynomial to use as the basis. This argument is optional,
%            by default, the monomials are used. The polynomials supported
%            by the vandermonde function are supported here.
%   origin - Specify an explicit value for the function at the origin. This
%            argument is optional, by default the intercept parameter is
%            left free for the approximation to choose.
%   alg    - Arguments to pass to the MATLAB linprog solver through
%            optimoptions. This argument must be a cell array of the option
%            pairs.
%
% Outputs:
%   coeff - The coefficients for the polynomial
%   E     - The maximum approximation error over the points
%
%
% Created by: Ian McInerney
% Created on: January 30, 2018
% Version: 1.0
% Last Modified: January 30, 2018
%
% Revision History
%   1.0 - Initial release

%% Parse the input
p = inputParser;
addOptional(p, 'poly', 'Monomial', @ischar);
addOptional(p, 'origin', NaN);
addOptional(p, 'alg', {});
parse(p, varargin{:});

poly = p.Results.poly;
origin = p.Results.origin;
alg = p.Results.alg;

[nx, ~] = size(x);


%% Create the Vandermonde function to call
if ( strcmp(poly, 'SSChebyshev') == 1)
    van = @(x) vandermonde(x, n+1, poly, min(x(:,1)), max(x(:,1)));
else
    van = @(x) vandermonde(x, n+1, poly);
end


%% Create the Vandermonde matrix for the constraint
V = van( x(:,1) );

V1 = [V, -ones(nx, 1)];
V2 = [V, ones(nx, 1)];

% Take into account a fixed origin if needed
if ( isnan(origin) )
    % Let the constant parameter be free
    b = x(:,2);
    coeff = [];
else
    % Force a value at the origin of the polynomial
    % This is equivalent to fixing the constant parameter and removing
    % it from the other side
    V1 = V1(:, 2:end);
    V2 = V2(:, 2:end);
    b = x(:,2) - origin.*ones(nx, 1);
    coeff = 1;
end

% Figure out how many coefficients there are
[~, nf] = size(V1);

% Create the inequality constraint matrices
A = [ V1;
     -V2];
b = [ b;
     -b];

% Want to minimize the last variable (error)
f = zeros(nf, 1);
f(end) = 1;


%% Solve the linear program
options = optimoptions('linprog', alg{:});      % Setup the options
c = linprog(f, A, b, [], [], [], [], options);


%% Extract the coefficients and the error
coeff = [coeff;
         c(1:end-1)];
E = c(end);

end

