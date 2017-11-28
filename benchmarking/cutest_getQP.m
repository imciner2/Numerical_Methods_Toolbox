function [ QP ] = cutest_getQP( varargin )
%CUTEST_GETQP Get the QP representation of a CUTEst problem
%
% This function will return the QP representation of the CUTEst problem.
% The optimization problem returned is of the form
%
%   min x'Qx + c'x
%   s.t. 
%       Aeq*x= beq
%       Ale*x <= ble
%
%       lb <= x <= ub
%
% Usage:
%   [Q, c, Aeq, beq, Ale, ble, ub, lb, x0, mu0] = cutest_getQP()
%   [Q, c, Aeq, beq, Ale, ble, ub, lb, x0, mu0] = cutest_getQP( problemDir )
%
% Inputs:
%   problemDir - Directory the mcutest.mex file is located in (if not in
%                current directory)
%
% Outputs:
%   QP  - A data structure containing the problem matrices (as follows)
%         Q    - The quadratic term matrix
%         c    - The linear term vector
%         Aeq  - The coefficient matrix for the linear equality constraints
%         beq  - The constant vector for the linear equality constraints
%         Ale  - The coefficient matrix for the linear inequality constraints
%         ble  - The constant vector for the linear inequality constraints
%         ub   - The upper bounds on the variables
%         lb   - The lower bounds on the variables
%         x0   - Initial guess for x
%         mu0  - Initial geuss for mu
%         name - Problem name
%
%
% Created by: Ian McInerney
% Created on: November 15, 2017
% Version: 1.2
% Last Modified: November 28, 2017
%
% Revision History:
%   1.0 - Initial Release
%   1.1 - Added initial variable output
%   1.2 - Modified output to be in a data structure, added problem name field


%% Make sure the cutest problem exists on the path
if ( ~isempty(varargin) )
    % Navigate to the proper directory
    originalDir = cd( varargin{1} );
else
    originalDir = pwd;
end

if ( exist('mcutest') ~= 3 )
    error('Cutest problem is not in the current directory');
end

% Get the CUTEst problem structure from the mex file
Problem = cutest_setup();


%% Get the objective function pieces
point = zeros(Problem.n, 1);        % Create the point x=0
QP.Q = cutest_isphess(point, 0);       % Get the hessian at x=0

QP.c = sparse( cutest_grad(point) );   % Get the gradient at x=0


%% Get the constraints
[con_rhs, con_lhs ] = cutest_scons(point);

% Check for non-linear constraints and warn the user
if ( sum(Problem.linear) ~= Problem.m )
    warning('Problem contains non-linear constraints. Those constraints will be ignored.');
end

% Create logical indexing for the different constraint types
con_eq = Problem.linear & Problem.equatn;
con_ineq = Problem.linear & ~Problem.equatn;

% Pull out the equality constraints from the constraint arrays
QP.Aeq = con_lhs( con_eq, : );
QP.beq = con_rhs( con_eq );

% Find the constraints that are LEQ or GEQ
% When the bound doesn't exist, it is put to +/-1e20 in the array
con_leq = con_ineq & (Problem.cl == -1e20);
con_geq = con_ineq & (Problem.cu == 1e20);

% Pull out the <= constraints
Ale = con_lhs( con_leq, : );
ble = Problem.cu( con_leq );

% Pull out the >= constraints and convert them to <=, then append them to
% the constraint matrix
QP.Ale = [Ale;
       -1.*con_lhs( con_geq, : ) ];
QP.ble = [ble;
       -1.*Problem.cl( con_geq) ];


%% Get the upper and lower bounds for the variables
QP.ub = Problem.bu;
QP.lb = Problem.bl;


%% Get the initial guesses for the variables
QP.x0 = Problem.x;
QP.mu0 = Problem.v;


%% Terminate the CUTEst problem before exiting
QP.name = Problem.name;
cutest_terminate();

cd(originalDir);

end
