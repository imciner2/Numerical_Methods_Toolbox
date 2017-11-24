function [D1, D2, numIter] = equil_simul_ruiz( A, varargin )
%EQUIL_SIMUL_RUIZ Perform simultaneous equilibration of rows and columns
%
% This function will perform simultaneous equilibration of the rows and
% columns of the matrix A. It is based on algorithm 2.1 given by Ruiz
% inside the report RAL-TR-2001-034
% (http://www.numerical.rl.ac.uk/reports/drRAL2001034.pdf)
%
% This algorithm will equilibrate the matrix in the infinity norm (so
% the infinity norm of each row and column is 1).
%
% The stopping criteria for this algorithm is to stop when the maximum
% value of every row and column is within a certain tolerance of 1. By
% default, this tolerance is 1e-4, but it can be changed.
%
% Usage:
%   [D1, D2] = equil_simul_ruiz( A )
%
% Inputs:
%   A - The matrix to perform equilibration on
%   Optional:
%       'Tolerance' - Specify the stopping tolerance
%
% Outputs:
%   D1 - The left-multiply matrix (scaling of the rows)
%   D2 - The right-multiply matrix (scaling of the columns)
%
%
% Created by: Ian McInerney
% Created on: November 21, 2017
% Version: 1.0
% Last Modified: November 21, 2017
%
% Revision History:
%   1.0 - Initial Release


%% Parse some input
if ( mod(length(varargin), 2) ~= 0)
    error('Incorrect optional arguments');
end

TOL = 1e-4;
if ( ~isempty(varargin) )
    for i=1:2:length(varargin)
        switch( varargin{i} )
        case 'Tolerance'
            TOL = varargin{i+1};
        end
    end
end


%% Create the sparse left/right matrices
[m, n] = size(A);
if ( issparse(A) )
    D1 = speye(m);
    D2 = speye(n);
else
    D1 = eye(m);
    D2 = eye(n);
end

%% Create the vectors to reference the diagonals of the matrices
DcIndices = 1:1:m;
DrIndices = 1:1:n;

TOL = 1e-4;

STOP = 0;
numIter = 1;
while (~STOP)
    % Find the maximum in each row/column
    rmax = max(A, [], 2);
    cmax = max(A, [], 1);
    
    % Create the scaling matrices
    Dr = sparse(DrIndices, DrIndices, 1./sqrt( rmax ) );
    Dc = sparse(DcIndices, DcIndices, 1./sqrt( cmax ) );
    
    % Update A
    A = Dr*A*Dc;

    % Update the final scaling matrices
    D1 = D1*Dr;
    D2 = D2*Dc;
    
    % Check the termination condition
    if ( (max( 1 - rmax ) <= TOL) || (max( 1 - cmax ) <= TOL) )
        STOP = 1;
    end
    numIter = numIter + 1;
    
end

end