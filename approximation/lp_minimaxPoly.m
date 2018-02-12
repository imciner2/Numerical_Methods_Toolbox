function [ coeff, E, p ] = lp_minimaxPoly( n, x, varargin )
%LP_MINIMAXPOLY Compute the minimax polynomial using the primal LP
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
%   [ coeff, E, p ] = LP_MINIMAXPOLY( n, x );
%   [ coeff, E, p ] = LP_MINIMAXPOLY( n, x, poly );
%   [ coeff, E, p ] = LP_MINIMAXPOLY( n, x, poly, origin );
%   [ coeff, E, p ] = LP_MINIMAXPOLY( n, x, poly, origin, alg );
%
% Inputs:
%   n      - The order of the polynomial to fit with
%   x      - The points to use to do the fitting
%   poly   - The polynomial to use as the basis. This argument is optional,
%            by default, the monomials are used. The polynomials supported
%            by the vandermonde function are supported here.
%   origin - Specify an explicit value for the function at the origin. This
%            argument is optional, by default the intercept parameter is
%            left free for the approximation to choose. To not specify an
%            intercept, pass in NaN.
%   alg    - Arguments to pass to the MATLAB linprog solver through
%            optimoptions. This argument must be a cell array of the option
%            pairs.
%
% Outputs:
%   coeff - The coefficients for the polynomial
%   E     - The maximum approximation error over the points
%   p     - The points from the input set where the error is largest
%
%
% see also LP_MINIMAXPOLY_DUAL
%
% Created by: Ian McInerney
% Created on: January 30, 2018
% Version: 1.1
% Last Modified: February 9, 2018
%
% Revision History
%   1.0 - Initial release
%   1.1 - Updated problem formulation and error handling

%% Parse the input
ip = inputParser;
addOptional(ip, 'poly', 'Monomial', @ischar);
addOptional(ip, 'origin', NaN);
addOptional(ip, 'alg', {});
parse(ip, varargin{:});

poly = ip.Results.poly;
origin = ip.Results.origin;
alg = ip.Results.alg;

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
    coeff = [1];
    O = zeros(1, n);
end

% Create the inequality constraint matrices (note, these are for Ax>=b)
A_pi = [ V, e;
        -V, e;
         O, 1];
b_pi = [ b;
        -b;
         0];

 % Figure out how many coefficients there are
[~, nf] = size(A_pi);
 
% Want to minimize the last variable (error)
f = zeros(nf, 1);
f(end) = 1;


%% Solve the linear program
options = optimoptions('linprog', alg{:});      % Setup the options
[c, pval, evalue] = linprog(f, -A_pi, -b_pi, [], [], [], [], options);


%% Extract the coefficients and the error
if (evalue == 1)
    coeff = [coeff;
             c(1:end-1)];
    E = c(end);
else
    warning('lp_minimaxPoly:Failed to find minimax polynomial');
    E = NaN;
end


%% Extract the points where the largest error occurs
% (e.g. where the constraints are active)
cval = abs(A_pi*c - b_pi);
cval = [ [1:1:nx, 1:1:nx]', cval(1:end-1)];
cval = sortrows(cval, 2);

ind = cval(1:n+1, 1);
p = x(ind);
p = sort(p);

end

