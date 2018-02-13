function [ A ] = eigTestMat( n, rho, lam1, lamn, varargin )
%EIGTESTMAT Create a square test matrix with a specified distribution of
%eigenvalues
%
% This function will generate a random matrix of size n where the
% eigenvalues are distributed according to equation 4 of 
%   A. Greenbaum and Z. Strakos, “Predicting the Behavior of Finite
%   Precision Lanczos and Conjugate Gradient Computations,” SIAM Journal
%   of Matrix Analysis an Applications, vol. 13, no. 1, pp. 121–137,
%   Jan. 1992.
%
% This distribution is modified by the parameter rho, which can vary
% between rho=0 (when the eigenvalues will be clustered at lam1 with an
% outlier at lamn), and rho=1 (where the eigenvalues will be equispaced
% between lam1 and lamn).
%
%
% Usage:
%   [ A ] = EIGTESTMAT( n, rho, lam1, lamn );
%   [ A ] = EIGTESTMAT( n, rho, lam1, lamn, U );
%
% Inputs:
%   n    - The size of the matrix
%   rho  - The distribution parameter
%   lam1 - The smallest eigenvalue
%   lamn - The largest eigenvalue
%   U    - The orthogonal eigenvectors (optional argument, if not specified
%          then a random orthogonal matrix will be generated).
%
% Outputs:
%   A - The matrix
%
%
% Created by: Ian McInerney
% Created on: February 13, 2018
% Version: 1.0
% Last Modified: February 13, 2018
%
% Revision History
%   1.0 - Initial release



%% Create the random eigenvectors
[U, ~] = qr(rand(n));


%% See if the eigenvectors were passed as input
p = inputParser;
addOptional(p, 'U', U);

parse(p, varargin{:});

U = p.Results.U;


%% Create the eigenvalues
lam = ones(n, 1);
lam(1) = lam1;
for i=2:1:n-1
    lam(i) = lam1 + (i-1)/(n-1)*(lamn - lam1)*rho^(n-i);
end
lam(n) = lamn;


%% Create the matrix
A = U*diag(lam)*U';

end

