function [f_all,gnorm_all] = gradmeth(fun, x0, tol, maxit)
% Gradient method with backtracking line search.
%
% Authour: Terrence Alsup
% Date: March 10, 2020
% File: gradmeth.m

x = x0;
[f, g] = fun(x); % Get the starting values.
gnorm = norm(g); 

f_all(1) = f;
gnorm_all(1) = gnorm;

it = 1;
% Check the stopping condition after each iteration.
while it <= maxit && gnorm > tol
    t = linesearch(fun, x, -g); % Backtracking linesearch for step length.
    x = x - t*g;                % The gradient descent update.
    
    [f, g] = fun(x);
    gnorm = norm(g);
    f_all(it + 1) = f;
    gnorm_all(it + 1) = gnorm;
    
    it = it + 1;
end

end