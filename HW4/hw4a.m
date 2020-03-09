% Solve the two-way partition SDP relaxation for the two datasets.
%
% Author: Terrence Alsup
% Date: March 4, 2020

% First solve for data set 1.

% Load the data.
W1 = open('hw4data1.mat').W;
W2 = open('hw4data2.mat').W;


n1 = size(W1, 1);
n2 = size(W2, 1);

% Solve for the first dataset.
cvx_begin sdp
    variable v(n1)
    maximize -sum(v)
    W1 + diag(v) >= 0;
cvx_end

d1 = cvx_optval

% Solve for the second dataset.
cvx_begin sdp
    variable v(n2)
    maximize -sum(v)
    W2 + diag(v) >= 0;
cvx_end

d2 = cvx_optval


% Solve the 2nd SDP for the first dataset.
cvx_begin sdp
    variable X(n1,n1) symmetric
    minimize trace(W1*X)
    diag(X) == 1;
    X >= 0;
cvx_end

d1_2 = cvx_optval