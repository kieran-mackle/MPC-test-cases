function output = forward_sim(A, B, C, D, xk, u)

% Currently for single output only - will need to change u indexing

n = size(A, 1);     % Number of states
m = size(B, 2);     % Number of control inputs
p = size(C, 1);     % Number of outputs
N = length(u);

x = zeros(n, N);
y = zeros(p, N);

x_0 = xk;
for i = 1:N
    y(:,i) = C*x_0 + D*u(i);
    x_new = A*x_0 + B*u(i);
    x_0 = x_new;
end



output = y;