function t = linesearch(fun, x, dx)
% Backtracking line search.

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