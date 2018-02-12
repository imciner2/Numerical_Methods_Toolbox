function [ coeff, E ] = lp_minimaxPoly_dual( n, x, varargin )
%LP_MINIMAXPOLY_DUAL Compute the minimax polynomial using the dual LP
%formulation
%
% Compute the minimax polynomial approximation of order n to a discrete set
% of points. This method uses the linear programming formulation of the
% problem, and solves it using MATLAB's linprog function (requires the
% Optimization Toolbox).
%
% The formulation for this problem is based upon the formulation presented
% in:
%   M. R. Osborne and G. A. Watson, “On the best linear Chebyshev
%   approximation,” The Computer Journal, vol. 10, no. 2, pp. 172–177, 1967.
%
%
% Usage:
%   [ coeff, E ] = LP_MINIMAXPOLY_DUAL( n, x );
%   [ coeff, E ] = LP_MINIMAXPOLY_DUAL( n, x, poly );
%   [ coeff, E ] = LP_MINIMAXPOLY_DUAL( n, x, poly, origin );
%   [ coeff, E ] = LP_MINIMAXPOLY_DUAL( n, x, poly, origin, alg );
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
% see also LP_MINIMAXPOLY
%
% Created by: Ian McInerney
% Created on: January 30, 2018
% Version: 1.0
% Last Modified: February 9, 2018
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
if ( strncmp(poly, 'SS', 2) == 1)
    % If the polynomial is a shifted and scaled version
    van = @(x) vandermonde(x, n+1, poly, min(x(:,1)), max(x(:,1)));
else
    % Normal polynomial
    van = @(x) vandermonde(x, n+1, poly);
end


%% Create the Vandermonde matrix for the constraint
V = van( x(:,1) );
e = ones(nx, 1);

% Take into account a fixed origin if needed
if ( isnan(origin) )
    % Let the constant parameter be free
    b = x(:,2);
    coeff = [];
    O = zeros(1, n+1);
else
    % Force a value at the origin of the polynomial
    % This is equivalent to fixing the constant parameter and removing
    % it from the other side
    V = V(:, 2:end);
    b = x(:,2) - origin.*ones(nx, 1);
    coeff = 1;
    O = zeros(1, n);
end

% Create the inequality constraint matrices
A_di = [e',  e'];
b_di = 1;

% Create the equality constraint matrices
A_de = [V', -V'];
b_de = O';

% Create the cost function
f_d = [b', -b'];


%% Solve the linear program
options = optimoptions('linprog', alg{:});      % Setup the options
[c, dval, evalue] = linprog(-f_d, A_di, b_di, A_de, b_de, zeros(length(f_d),1), [], options);


%% Extract the coefficients and the error
if (evalue == 1)
    coeff = [coeff;
             c(1:end-1)];
    E = dval;
else
    warning('lp_minimaxPoly_dual:Failed to find minimax polynomial');
    E = NaN;
end

end

