function [p, v, Z] = solve_sdp1(W)
% Solve the dual of the two-way partitioning problem and return the 
% optimal value p, solution v, and dual variable Z.
%
% Author: Terrence Alsup
% Date: March 4, 2020

n = size(W, 1); % Get the dimension.

cvx_begin sdp
variable v(n)
dual variable Z
maximize -sum(v)
W + diag(v) >= 0 : Z; % Z is the dual variable for this constraint.
cvx_end

p = cvx_optval;

end

