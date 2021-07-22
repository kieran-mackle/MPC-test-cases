% ======================================================================= %
%                           Constrained Ball MPC                          %
% ======================================================================= %

clearvars;

% ----------------------------------------------------------------------- %
% Define MPC Parameters
% ----------------------------------------------------------------------- %
Ts = 0.5;             % Sampling interval
Hp = 10;            % Prediction horizon
total_time = 50;    % Total time horizon

% ----------------------------------------------------------------------- %
% Define initial state
% ----------------------------------------------------------------------- %
x0  = [2, -0.5]';
u0  = [-0.5]';
z0  = ball_output(x0);

% ----------------------------------------------------------------------- %
% Define LTI model
% ----------------------------------------------------------------------- %
i = 1;              % Model index
% Construct input from initial state
input.phase.time    = 0;
input.phase.state   = x0;
input.phase.control = u0;


% Call myJac.m to evaluate the Jacobians at the input time
A = myJac(@ball_dynamics, input, 'state');
B = myJac(@ball_dynamics, input, 'control');

% Construct C and D
C = numerical_jac(@ball_output, input.phase.state);
D = zeros(size(C,1), size(B,2));

n = size(A, 1);     % Number of states
m = size(B, 2);     % Number of control inputs
p = size(C, 1);     % Number of outputs

% ----------------------------------------------------------------------- %
% Simulate trajectory of ball
% ----------------------------------------------------------------------- %
N           = round(total_time/Ts);
t           = zeros(1, N+1);
U           = zeros(m, N+1);
X           = zeros(n, N+1);
Z           = zeros(p, N+1);
set_points  = zeros(p, N+1);
r           = [5]';

% Fill initial state, control and output
X(:,1)      = x0;
U(:,1)      = u0;
Z(:,1)      = z0;

% % Initialise figures
% figure(1);
% clf;
% title('Control inputs');
% 
% figure(2); 
% clf;
% sgtitle('MPC predicted trajectory');
% subplot(2,1,1);
% hold on;
% grid on;
% title('Position');
% ylabel('x (m)');
% subplot(2,1,2);
% hold on;
% grid on;
% title('Velocity');
% ylabel('v (m/s)');
% xlabel('t (s)');


for k = 1:N
    % Define plant and convert to discrete-time domain
    % ------------------------------------------------
    plant   = ss(A, B, C, D);
    plant   = c2d(plant, Ts);
    Ad      = plant.A;
    Bd      = plant.B;
    
    plant2  = ss(A, eye(size(A)), C, D);
    plant2  = c2d(plant2, Ts);
    
    % Define current state and transform system
    % -----------------------------------------
    x0          = input.phase.state;       % x(k)
    u0          = input.phase.control;     % u(k-1)
    xbar_k      = zeros(size(x0));
    ubar_km1    = zeros(size(u0));

    % Calculate current operating point
    % -----------------------------------------
    f0      = ball_dynamics(input).dynamics';
    f0d     = plant2.B * f0;
    g0      = ball_output(x0)';
    
    % Define cost weighting matrices
    % -----------------------------------------
    Q       = eye(p);                       % Output weightings
    R       = eye(m);                       % Control weightings

    % Define constraint matrices
    % -----------------------------------------
    % Input rate constraints (independent)
    E       = [-1,   -1;
                1,   -1];
    
    % Input constraints (relative to ubar)
    F       = [-1, -3 - u0(1);
               1,  -3 + u0(1)];

    % Output constraints (on physical system)
    G_const = [-1,  -10;
               1,  -10];

    
    % Format constraint and weighting matrices
    % -----------------------------------------
    Q_cell  = repmat({Q}, 1, Hp);           % Repeats Q Hp+1 times
    Q_fancy = blkdiag(Q_cell{:});           % Diagonal matrix of Qc
    R_cell  = repmat({R}, 1, Hp);
    R_fancy = blkdiag(R_cell{:});
    
    a       = size(E,1);
    E_nom   = E(1:a, 1:m);
    w_nom   = E(:, end);
    EE      = zeros(Hp*a, Hp*m+1);
           
    b       = size(F,1);
    F_nom   = F(1:b, 1:m);
    f_nom   = F(:, end);
    FF      = zeros(Hp*b, Hp*m+1);
    F_fancy = zeros(Hp*b, m*Hp);
    
    c       = size(G_const,1);
    g_nom   = G_const(:, end);
    G_nom   = G_const(1:c, 1:p);
    GG      = zeros(Hp*c, Hp*p+1);
    
    % Build prediction matrices
    % -----------------------------------------
    for i = 1:Hp
        EE(a*(i-1)+1:a*i, m*(i-1)+1:m*i)    = E_nom;
        EE(a*(i-1)+1:a*i, Hp*m+1)           = w_nom;

        FF(b*(i-1)+1:b*i, m*(i-1)+1:m*i)    = F_nom;
        FF(b*(i-1)+1:b*i, Hp*m+1)           = f_nom;

        GG(c*(i-1)+1:c*i, p*(i-1)+1:p*i)    = G_nom;
        GG(c*(i-1)+1:c*i, Hp*p+1)           = g_nom;
    end

    % ------------------------------------------------------------------- %
    % Build Prediction Matrices
    % ------------------------------------------------------------------- %
    AB_big      = zeros(n*Hp, m);
    ABdu_big    = zeros(n*Hp, m*Hp);
    T_fancy     = zeros(p*Hp, 1);
    IA_big      = zeros(n*Hp, n);

    for i = 1:Hp
        A_big(n*(i-1)+1:n*i, 1:n)           = Ad^i;
        C_big(p*(i-1)+1:p*i, n*(i-1)+1:n*i) = C;

        AB_sum = 0;
        for j = 0:i-1
            AB_sum = AB_sum + Ad^j * Bd;
        end
        AB_big(n*(i-1)+1:i*n, 1:m) = AB_sum;

        ABdu_big(n*(i-1)+1:n*i, 1:m) = AB_sum;
        if i > 1
            ABdu_big(n*(i-1)+1:n*i, m+1:m*Hp) = ABdu_big(n*(i-2)+1:n*(i-1), 1:m*(Hp-1));
        end

        T_fancy(p*(i-1)+1:p*i, 1) = r;

        % Constraints
        F_sum = zeros(b*Hp, m);
        for j = i:Hp
            F_sum = F_sum + FF(:, m*(j-1)+1:m*j);
        end
        F_fancy(:, m*(i-1)+1:m*i) = F_sum;
        
        IA_sum = 0;
        for j = 0:i-1
            IA_sum = IA_sum + Ad^j;
        end
        IA_big(n*(i-1)+1:i*n, 1:n) = IA_sum;

    end

    W       = EE(:, 1:end-1);
    w       = EE(:, end);
    f       = FF(:,end);
    Tau     = GG(:, 1:end-1);
    g       = GG(:, end);
    
    F_1     = F_fancy(:, 1:m);
    psi     = C_big*A_big;
    gamma   = C_big*AB_big;
    theta   = C_big*ABdu_big;
    phi     = C_big*IA_big;
    
    % Formulate Quadratic Program
    % ----------------------------------
    % Prediction matrices
    epsilon = T_fancy - psi*xbar_k - gamma*ubar_km1 - phi*f0d - repmat(g0, Hp, 1);
    G       = 2 * theta' * Q_fancy * epsilon;
    H       = theta' * Q_fancy * theta + R_fancy;
    
    % Constraint matrices
    Omega   = [F_fancy; 
               Tau*theta; 
               W];
    omega   = [-F_1*ubar_km1 - f;
               -Tau*(psi*xbar_k + gamma*ubar_km1 + phi*f0d + repmat(g0, Hp, 1)) - g;
               -w];
           
    % Solve Quadratic Program
    % ----------------------------------
    [dUbar, fval, ~, ~, ~] = quadprog(2*H, -G, Omega, omega);
    
    % Get optimal control input
    % ----------------------------------
    Ubar                    = zeros(size(dUbar));
    ubar_previous           = ubar_km1;
    for h = 1:Hp
        delta_ubar              = dUbar(m*(h-1)+1:h*m,:);
        Ubar(m*(h-1)+1:h*m, 1)  = ubar_previous + delta_ubar;
        ubar_previous           = Ubar(m*(h-1)+1:h*m, 1);   
    end
    U_k     = Ubar + repmat(u0, Hp,1);
    uk      = U_k(1:m,:);
    
    % Calculate next state
    % ----------------------------------
    input.phase.time        = t(k);
    xk1                     = x0;
    for j = 1:100
        input.phase.state   = xk1';
        input.phase.control = uk';
        dynamics            = ball_dynamics(input);
        f_dot               = dynamics.dynamics;
        xk1                 = xk1 + f_dot'*Ts/100;
    end
    zk          = ball_output(xk1);
    
    % Calculate MPC predicted output
    % ----------------------------------
    z_predicted     = psi*xbar_k + gamma*ubar_km1 + theta*dUbar + ...
                        phi*f0d + repmat(g0, Hp, 1);
    xbar_predicted  = A_big*xbar_k + AB_big*ubar_km1 + ABdu_big*dUbar + ...
                        IA_big*f0d;
    x_predicted     = xbar_predicted + repmat(x0, Hp, 1);
    
    % Append results to output arrays
    % ----------------------------------
    t(k+1)              = Ts*k;
    U(:, k+1)           = uk;
    X(:, k+1)           = xk1;
    Z(:, k+1)           = zk;
    set_points(:, k+1)  = r;
    
    
    % Plotting
    % ----------------------------------
%     figure(1);
%     hold on; grid on;
%     stairs(linspace(0, Hp*Ts, Hp+1), [u0(1); U_k(1:m:end)], 'k-');
%     ylabel('Input force (N)');
%     
%     figure(2);
%     subplot(2,1,1);
%     grid on;
%     plot(linspace(0, Hp*Ts, Hp+1), [z0(1); z_predicted(1:p:end)], 'k-');
%     subplot(2,1,2);
%     plot(linspace(0, Hp*Ts, Hp+1), [x0(2); x_predicted(2:n:end)], 'k-');
    
    % Re-calculate Jacobians
    % ----------------------------------
    input.phase.time    = t(k+1);
    input.phase.state   = xk1;
    input.phase.control = uk;
    
    A = myJac(@ball_dynamics, input, 'state');
    B = myJac(@ball_dynamics, input, 'control');
    C = numerical_jac(@ball_output, input.phase.state);
    D = zeros(size(C,1), size(B,2));
    
end


% ----------------------------------------------------------------------- %
% Extract results
% ----------------------------------------------------------------------- %
figure(1);
clf;
hold on;
grid on;
title('MPC input');
stairs(t, U(1,:))

figure(2);
clf;
sgtitle('Physical state trajectory');
subplot(2,1,1);
hold on;
grid on;
plot(t, X(1:n:end), 'k-');
ylabel('x-position (m)');
subplot(2,1,2);
hold on;
grid on;
plot(t, X(2:n:end), 'k-');
ylabel('velocity (m)');
xlabel('Time (s)');

figure(3);
clf;
hold on;
title('Output');
grid on;
plot(t, Z);
ylabel('x-position (m)');
xlabel('Time (s)');

