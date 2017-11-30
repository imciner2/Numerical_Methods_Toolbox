classdef AugmentedLagrangian_Hessian
    
    properties
        A
        C
        rho
    end
    methods
        function obj = AugmentedLagrangian_Hessian(A, C, rho)
            obj.A = A;
            obj.C = C;
            obj.rho = rho;
        end
        
        function res = mtimes(a, b)
            leftType = isa(a, 'AugmentedLagrangian_Hessian');
            rightType = isa(b, 'AugmentedLagrangian_Hessian');
            
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