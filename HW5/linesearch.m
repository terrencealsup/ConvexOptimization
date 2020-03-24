function t = linesearch(fun, x, dx)
% Backtracking line search.
% Given a function fun, a point x, and a direction dx, compute the 
% step length to take.
%
% Author: Terrence Alsup
% Date: March 10, 2020
% File: linesearch.m

% Fixed parameters for the backtracking line search.
a = 0.25;
b = 0.5;

t = 1; % Define the initial step length.
[fx, gx] = fun(x); % Evaluate the function at the initial point.

% Do the backtracking.
while fun(x + t*dx) > fx + a*t*gx'*dx
    t = t * b;
end

end