function [ Q, c, Aeq, beq, Ale, ble, ub, lb ] = cutest_getQP( )
%CUTEST_GETQP Get the QP representation of a CUTEst problem
%
% This function will return the QP representation of the CUTEst problem.
% The optimization problem returned is of the form
%
%   min x'Qx + c'x
%   s.t. 
%       Aeq*x= b
%       Ale*x <= ble
%
%       lb <= x <= ub
%
% Outputs:
%   Q   - The quadratic term matrix
%   c   - The linear term vector
%   Aeq - The coefficient matrix for the linear equality constraints
%   beq - The constant vector for the linear equality constraints
%   Ale - The coefficient matrix for the linear inequality constraints
%   ble - The constant vector for the linear inequality constraints
%   ub  - The upper bounds on the variables
%   lb  - The lower bounds on the variables
%
%
% Created by: Ian McInerney
% Created on: November 15, 2017
% Version: 1.0
% Last Modified: November 16, 2017
%
% Revision History:
%   1.0 - Initial Release


%% Make sure the cutest problem exists on the path
if ( exist('mcutest') ~= 3 )
    error('Cutest problem is not loaded into path. Make sure the mcutest mex file is on the path.');
end

% Get the CUTEst problem structure from the mex file
Problem = cutest_setup();


%% Get the objective function pieces
point = zeros(Problem.n, 1);        % Create the point x=0
Q = cutest_isphess(point, 0);       % Get the hessian at x=0

c = sparse( cutest_grad(point) );   % Get the gradient at x=0


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
Aeq = con_lhs( con_eq, : );
beq = con_rhs( con_eq );

% Find the constraints that are LEQ or GEQ
% When the bound doesn't exist, it is put to +/-1e20 in the array
con_leq = con_ineq & (Problem.cl == -1e20);
con_geq = con_ineq & (Problem.cu == 1e20);

% Pull out the <= constraints
Ale = con_lhs( con_leq, : );
ble = Problem.cu( con_leq );

% Pull out the >= constraints and convert them to <=, then append them to
% the constraint matrix
Ale = [Ale;
       -1.*con_lhs( con_geq, : ) ];
ble = [ble;
       -1.*Problem.cl( con_geq) ];


%% Get the upper and lower bounds for the variables
ub = Problem.bu;
lb = Problem.bl;


%% Terminate the CUTEst problem before exiting
cutest_terminate();

end
