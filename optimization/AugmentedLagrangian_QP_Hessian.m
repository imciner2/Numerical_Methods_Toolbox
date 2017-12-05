classdef AugmentedLagrangian_QP_Hessian
    %AugmentedLagrangian_QP_Hessian The Hessian matrix for an Augmented Lagrangian QP problem
    %
    % The Augmented Lagrangian for QPs solves the problem of the form:
    %         min 1/2 x'Ax - b'x
    %         s.t. Cx = d
    % by introducing a quadratic penalty term and also the Lagrange
    % multipler vector for the constraint. This then creates a saddle-point
    % optimization problem:
    %         max min 1/2 x'Ax - b'x u'*(Cx - d) + rho*(Cx - d)^2
    %          u   x
    % The Hessian of this optimization problem (which this class
    % represents) is:
    %         A + rho*C'C
    %
    % Properties:
    %   A   - The Hessian of the objective function
    %   C   - The LHS matrix of the constraints
    %   rho - The penalty parameter
    %
    % Methods:
    %   primalFunctional - Compute the Hessian of the primal functional
    %   value - Compute the value of the Hessian
    properties
        A       % The Hessian of the objective function
        C       % The LHS matrix of the constraints
        rho     % The penalty parameter
    end
    methods
        function obj = AugmentedLagrangian_QP_Hessian(A, C, rho)
            %AugmentedLagrangian_QP_Hessian Create the Hessian
            %
            %  Creates the Hessian object.
            %
            % Usage:
            %  [ Hess ] = AugmentedLagrangian_QP_Hessian(A, C, rho)
            %
            % Inputs:
            %  A   - The Hessian of the objective function
            %  C   - The LHS matrix of the constraints
            %  rho - The penalty parameter
            %
            % Outpus:
            %  Hess - The Hessian object
            obj.A = A;
            obj.C = C;
            obj.rho = rho;
        end
        
        function res = primalFunctional(a)
            %primalFunctional Compute the Hessian of the primal functional
            %
            %  Compute the Hessian of the primal functional using the
            %  formula given by Bertsekas in "Constrained Optimization
            %  and Lagrange Multiplier Methods", 1982. Namely
            %      hess(p) = (C (A)^(-1) C')^(-1)
            res = inv( (a.C) * inv((a.A)) * (a.C)');
        end
        
        function res = value(a)
            %value Evaluate the Hessian matrix and return the matrix
            
            res = a.A + a.rho.*a.C'*a.C;
        end
        
        function res = subsref(a, b)
            %subsref Get the Hessian value at specific indices
            
            % Test the type of reference to make sure it is right
            if ( ~strcmp(b.type, '()') )
                error('Unsupported reference type');
            end
            
            % Make sure there are enough subscripts
            [t1, t2] = size(b.subs);
            if ( (t1 ~= 1) && (t2 ~= 2) )
                error('Not enough subscripts');
            end
            
            % Get the relevant indices to pull from
            [Arow, Acol] = size(a.A);
            if ( b.subs{1} == ':' )
                % It is all the elements in the column
                rowIndices = 1:1:Arow;
            else
                % It is only specific elements in the column
                rowIndices = b.subs{1};
            end
            
            if ( b.subs{2} == ':' )
                % It is all the elements in the column
                colIndices = 1:1:Acol;
            else
                % It is only specific elements in the column
                colIndices = b.subs{2};
            end
            
            % Create the transpose of C
            cTrans = (a.C)';
            
            % Actually create the result 
            res = a.A(rowIndices, colIndices);
            res = res + (a.rho).*cTrans(rowIndices, :)*(a.C(:, colIndices));

        end
        
        function res = transpose(a)
            %transpose Compute the transpose of the Hessian
            res = AugmentedLagrangian_QP_Hessian( (a.A).', (a.C), a.rho);
        end
        
        function res = ctranspose(a)
            %ctranspose Compute the conjugate transpose of the Hessian
            res = AugmentedLagrangian_QP_Hessian( (a.A)', (a.C), a.rho);
        end
        
        function res = mtimes(a, b)
            %mtimes Compute the matrix multiplication of the Hessian
            leftType = isa(a, 'AugmentedLagrangian_QP_Hessian');
            rightType = isa(b, 'AugmentedLagrangian_QP_Hessian');
            
            if ( leftType && ~rightType )
                % Left multiply a vector by the Hessian
                % Compute (A + rho*C'*C)*b with b a vector
                
                % Compute A*b
                A_res = (a.A)*b;
                
                % Compute (rho*C'*C)*b from right-to-left
                C_res = (a.C)*b;
                C_res = (a.C)'*C_res;
                C_res = (a.rho).*C_res;     % rho is a scalar
                
                % Compute the final result
                res = A_res + C_res;
                
            elseif ( ~leftType && rightType )
                % Right multiply a vector by the Hessian
                % Compute a*(A + rho*C'*C) with a a vector
                
                % Compute a*A
                A_res = a*(b.A);
                
                % Compute a*(rho*C'*C) from left-to-right
                C_res = a.*(b.rho);         % rho is a scalar
                C_res = C_res*(b.C)';
                C_res = C_res*(b.C);
                
                % Compute the final result
                res = A_res + C_res;
                
            elseif ( leftType && rightType )
                % Multiply a Hessian with a Hessian
                
            else
                % This is bad...
                res = NaN;
            end
        end
        
    end
end