function [mincut, maxcut, x] = partition(W)
% Brute force computation of the two-way partition problem.
% mincut is the optimal value for x'Wx
% maxcut is the optimal max cut value from Goemans and Williamson.
% x is the vector of +-1 according to the partition.
%
% Author: Terrence Alsup
% Date: March 4, 2020

n = size(W, 1); % Get the dimension.

% Loop over all 2^n possible values.
x = ones(n, 1);
mincut = x'*W*x;
for i=0:2^n - 1
    % Compute the decimal number to a binary (row) vector.
    temp_x = de2bi(i, n);
    % Convert the 0's to -1's.
    temp_x(temp_x == 0) = -1;
    temp_cut = temp_x*W*temp_x';
    if temp_cut < mincut
        x = temp_x';
        mincut = temp_cut;
    end
end

% Now compute the max cut value from the Goemans and Williamson paper.
maxcut = 0.25 * (ones(n,1)'*W*ones(n,1) - mincut);

end

