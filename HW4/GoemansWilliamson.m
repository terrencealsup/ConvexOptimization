function [mincut, maxcut, x] = GoemansWilliamson(W)
% Partition assignment according to Goemans and Williamson.
% First solve the relaxed SDP and then do the assignment.
%
% mincut is the value x'Wx according to the partition.
% maxcut is the value (1'W1 - x'Wx)/4 from Goemans and Williamson.
% x is the partition.
%
% Author: Terrence Alsup
% Date: March 4, 2020

n = size(W, 1);

% Solve the SDP relaxation for X.
cvx_begin sdp
variable X(n,n) symmetric
dual variables Z y
minimize trace(W*X)
X >= 0 : Z;         % Z is a matrix and a dual variable.
diag(X) == 1 : y;   % y is a vector and the dual variable here.
cvx_end

% Compute V via Cholesky and make sure that it is positive definite.
V = chol(X + 1e-12*eye(n));

% The vector r on the unit sphere.
r = ones(n,1)/sqrt(n); 

% Do the assignment.
x = zeros(n, 1);
for i = 1:n
    if V(:,i)'*r >= 0
        x(i) = 1;
    else
        x(i) = -1;
    end
end

% Compute the cut-values.
mincut = x'*W*x;
maxcut = 0.25 * (ones(n,1)'*W*ones(n,1) - mincut);

end

