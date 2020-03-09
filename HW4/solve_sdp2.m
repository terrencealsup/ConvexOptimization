function [p, X, Z, y] = solve_sdp2(W)
% Solve the SDP relaxation of the two-way partitioning problem and return
% the optimal value p, solution X, and dual variables Z, y.
%
% Author: Terrence Alsup
% Date: March 4, 2020

n = size(W, 1); % Get the dimension.

cvx_begin sdp
variable X(n,n) symmetric
dual variables Z y
minimize trace(W*X)
X >= 0 : Z;         % Z is a matrix and a dual variable.
diag(X) == 1 : y;   % y is a vector and the dual variable here.
cvx_end

p = cvx_optval;
end

