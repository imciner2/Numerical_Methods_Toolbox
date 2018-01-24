% This script will test the remes_exchange function using example 2 on
% page 301 of 
%   W. Fraser, “A Survey of Methods of Computing Minimax and Near-Minimax
%   Polynomial Approximations for Functions of a Single Independent
%   Variable,” J. ACM, vol. 12, no. 3, pp. 295–314, 1965.


%% Create the function value pairs
x  = 0:0.2:3;
fx = [0.00000;
      0.44721;
      0.63245;
      0.77460;
      0.89443;
      1.00000;
      1.09545;
      1.18322;
      1.26491;
      1.34164;
      1.41421;
      1.48324;
      1.54919;
      1.61245;
      1.67332;
      1.73205];


%% Call the function
% No polynomial specified, defaults to monomial
fprintf('Testing default polynomial\n');
[ b, E ] = remes_exchange(3, [x', fx])

% Monomial basis
fprintf('\nTesting Monomial\n');
[ b, E ] = remes_exchange(3, [x', fx], 'Monomial')

% Chebyshev basis
fprintf('\nTesting Chebyshev\n');
[ b, E ] = remes_exchange(3, [x', fx], 'Chebyshev')