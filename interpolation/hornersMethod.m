function P = hornersMethod(x, F, xbar)
% HORNERSMETHOD Use Horner's method to interpolate points
%
% Calculate the value of an interpolated function using Horner's method to
% evaluate the Newton form of the interpolating polynomial.
%
%
% Usage:
%   [ P ] = HORNERSMETHOD(x, F, xbar)
%   
% Inputs:
%   x    - The x points the interpolation goes through
%   F    - The divided difference table for the function
%   xbar - The x coordinates of the points to find the value of
%
% Outputs:
%   P - Interpolated value of the xbar points
%
%
% see also DIVIDEDDIFFTABLE
%
% Created by: Ian McInerney
% Created on: September 21, 2016
% Version: 1.0
% Last Modified: February 12, 2018
%
% Revision History
%   1.0 - Initial release


%% Get the coefficients for the Newton forward form
coeff = diag(F);

%% Get some lengths
len = length(xbar);
numInterPoints = length(x);

%% Initialize P
P = nan(len,1);

%% Iterate over each xbar
for (i = 1:1:len)
    % Set the initial value to the last coefficient
    P(i) = coeff(numInterPoints);
    for (j = (numInterPoints-1):-1:1)
        % Update the output value with the next point
        P(i) = P(i)*(xbar(i) - x(j)) + coeff(j);
    end
end

end