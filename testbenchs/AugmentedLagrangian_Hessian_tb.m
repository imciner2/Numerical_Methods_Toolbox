numTests = 100;
numVar = 10;
numCon = 100;

tol = 1e-6;

rightTestPass = 0;
leftTestPass = 0;
transTestPass = 0;
ctransTestPass = 0;
subsTestPass = 0;

for (i = 1:1:numTests)
    % Randomly generate the A and C matrices
    A = rand(numVar, numVar);
    C = rand(numCon, numVar);
    
    rho = rand(1,1);
    
    % Create the Hessian
    Hess = AugmentedLagrangian_QP_Hessian(A, C, rho);
    
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
    
    % Test conjugate transpose
    test = value(Hess');
    truth = (A+rho*(C'*C))';
    
    if ( all(all(abs(test - truth) > tol)) )
        warning('Failed conjugate transpose test');
    else
        ctransTestPass = ctransTestPass + 1;
    end
    
    % Test transpose
    test = value(Hess.');
    truth = (A+rho*(C'*C)).';
    
    if ( all(all(abs(test - truth) > tol)) )
        warning('Failed transpose test');
    else
        transTestPass = transTestPass + 1;
    end
    
    % Test subscript referencing
    numRefs = randi(numVar);
    row = randi(numVar, 1, numRefs);
    col = randi(numVar, 1, numRefs);
    test = Hess(row, col);
    truth = (A+rho*(C'*C));
    truth = truth(row, col);
    
    if ( all(all(abs(test - truth) > tol)) )
        warning('Failed subscript referencing test');
    else
        subsTestPass = subsTestPass + 1;
    end
    
end

disp(['Right multiplication test pass ', num2str(rightTestPass), ' out of ', num2str(numTests)]);
disp(['Left multiplication test pass ', num2str(leftTestPass), ' out of ', num2str(numTests)]);
disp(['Conjugate transpose test pass ', num2str(ctransTestPass), ' out of ', num2str(numTests)]);
disp(['Transpose test pass ', num2str(transTestPass), ' out of ', num2str(numTests)]);
disp(['Subscript reference test pass ', num2str(subsTestPass), ' out of ', num2str(numTests)]);