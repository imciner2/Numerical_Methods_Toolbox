numTests = 100;
numVar = 10;
numCon = 40;

tol = 1e-6;

rightTestPass = 0;
leftTestPass = 0;

for (i = 1:1:numTests)
    % Randomly generate the A and C matrices
    A = rand(numVar, numVar);
    C = rand(numCon, numVar);
    
    rho = rand(1,1);
    
    % Create the Hessian
    Hess = AugmentedLagrangian_Hessian(A, C, rho);
    
    vec = rand(numVar, 1);
    
    % Test right multiplication
    test = vec'*Hess;
    truth = vec'*(A+rho*(C'*C));
    
    if ( abs(test - truth) > tol )
        warning('Failed right multiplication test');
    else
        rightTestPass = rightTestPass + 1;
    end
    
    % Test left multiplication
    test = vec'*Hess;
    truth = vec'*(A+rho*(C'*C));
    
    if ( abs(test - truth) > tol )
        warning('Failed left multiplication test');
    else
        leftTestPass = leftTestPass + 1;
    end
    
end

disp(['Right multiplication test pass ', num2str(rightTestPass), ' out of ', num2str(numTests)]);
disp(['Left multiplication test pass ', num2str(leftTestPass), ' out of ', num2str(numTests)]);