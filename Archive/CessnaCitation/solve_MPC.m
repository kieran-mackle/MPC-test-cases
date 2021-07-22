function [opt_u, opt_y, t] = solve_MPC(A,B,C,D,Hp,Hu,Ts, tolerance)

% Define conversion constants
deg = pi/180;
rad = 180/pi;

% Determine sizes
n = size(A, 1);     % Number of states
m = size(B, 2);     % Number of control inputs
p = size(C, 1);     % Number of outputs

r = [0, 40, 0]';    % Set-point - include as input
Q = eye(p);         % State weightings
R = eye(m);         % Control weightings

% Define plant and convert to discrete-time domain
plant = ss(A, B, C, D);
plant = c2d(plant, Ts);

Ad = plant.A;
Bd = plant.B;
Cd = plant.C;
Dd = plant.D;

% Create full simulation vectors
Cc = repmat({Cd}, 1, Hp+1);
CC = blkdiag(Cc{:});
Dc = repmat({Dd}, 1, Hp+1);
DD = blkdiag(Dc{:});
rc = repmat({r'}, 1, Hp+1);
RR = cell2mat(rc)';
Qc = repmat({Q}, 1, Hp+1);
QQ = blkdiag(Qc{:});
big_Rc = repmat({R}, 1, Hp);
big_R = blkdiag(big_Rc{:});

% Simulate Trajectory
% N = 20;
% t = zeros(1,N);
% opt_x = zeros(n, N);
% opt_y = zeros(p, N);
% opt_u = zeros(m, N);

x_0 = zeros(n,1); % Initial state vector
u_0 = zeros(m,1); % Initial control vector


f = waitbar(0);
x_k = x_0;
u_k = u_0;

cvx_optval = 5*tolerance;
k = 1;

while cvx_optval > tolerance
    waitbar(tolerance/cvx_optval, f, sprintf('Convergence (%d%%)', tolerance/cvx_optval*100));
    
    % Add results to vectors
    opt_u(:,k) = u_k;
    opt_x(:,k) = x_k;
    opt_y(:,k) = Cd*x_k + Dd*u_k;
    
    % Run cvx
    cvx_begin quiet
        variable x(n, Hp+1)
        variable u(m*(Hp+1), 1)

        XX = reshape(x, [n*(Hp+1),1]);
        UU = u;
        YY = CC*XX + DD*UU;

        du = u(m+1:end) - u(1:end-m);
        
        minimize( (YY-RR)' * QQ * (YY-RR) + du'*big_R*du)
        subject to
            x(:, 1) == x_k;
            -15*deg <= u(:) <= 15*deg % elevator angle
            -30*deg*Ts <= du(:) <= 30*deg*Ts % Elevator slew rate
            -30*deg*Ts <= u(m+1:2*m) - u(1:m) <= 30*deg*Ts 
                % Initial elevator slew rate

            
        for i = 1:Hu
            subject to
                x(:, i+1) == Ad * x(:, i) + Bd * u(m*(i-1)+1:m*i);
                -20*deg<= YY(p*(i-1)+1,1) <= 20*deg  % pitch angle
        end
        
        for i = Hu+1:Hp
            subject to
                u(m*(i-1)+1:m*i) == u(m*(i-2)+1:m*(i-1));
                x(:, i+1) == Ad * x(:, i) + Bd * u(m*(i-1)+1:m*i);
                -20*deg<= YY(p*(i-1)+1,1) <= 20*deg  % pitch angle
        end
        
    cvx_end
    
    % Extract latest results
    x_k1 = x(:,1);
    u_k1 = u(1:m);
    
    % Update state
    x_k = Ad*x_k1 + Bd*u_k1;
    u_k = u_k1;
    
    t(k) = Ts*(k-1);
    cvx_optval
    
    k = k+1;
end
close(f)




% output = [opt_u, opt_y];