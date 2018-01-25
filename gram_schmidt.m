function [ a ] = gram_schmidt( a, Q )
%GRAM_SCHMIDT Make the vector a orthogonal to the columns of Q
%
% Makes the vector a orthogonal to the columns of Q using the modified
% Gram-Schmidt algorithm.
%
%
% Usage:
%   [ a ] = gram_schmidt( a, Q );
%
% Inputs:
%   a - The vector to make orthogonal
%   Q - The matrix containing the other vectors as columns
%
% Outputs:
%   a - The orthogonal vector
%
%
% Created by: Ian McInerney
% Created on: January 25, 2018
% Version: 1.0
% Last Modified: January 25, 2018
%
% Revision History:
%   1.0 - Initial Release

[~, n] = size(Q);

for (i=1:1:n)
    alpha = a'*Q(:,i);
    beta = Q(:,i)'*Q(:,i);
    
    a = a - (alpha)/(beta)*Q(:,i);
end


end

