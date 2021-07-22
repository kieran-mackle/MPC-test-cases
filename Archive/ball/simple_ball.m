% ======================================================================= %
%                           Constrained Ball MPC                          %
% ======================================================================= %

clearvars;
deg = pi/180;
rad = 180/pi;

% ----------------------------------------------------------------------- %
% Define MPC Parameters
% ----------------------------------------------------------------------- %
params.timestep     = 0.5;
params.horizon      = 10;
params.sim_time     = 50;

% ----------------------------------------------------------------------- %
% Define function handles
% ----------------------------------------------------------------------- %
functions.dynamics  = @ball_dynamics;
functions.output    = @ball_output;

% ----------------------------------------------------------------------- %
% Define initial state
% ----------------------------------------------------------------------- %
x0  = [2, 0];
u0  = [0.5]';

initial.state       = x0;
initial.control     = u0;

target_ref          = [5]';

% ----------------------------------------------------------------------- %
% Define cost and constraint matrices
% ----------------------------------------------------------------------- %
cost.output         = eye(1);
cost.control        = eye(1);

constraints.hard.rate   = [-1, 1];
constraints.hard.input  = [-5, 5];
constraints.hard.output = [-10, 10];

constraints.soft.rate   = [-10, 10];
constraints.soft.input  = [-30, 30];
constraints.soft.output = [-10, 10];

% ----------------------------------------------------------------------- %
% Solve MPC QP Problem
% ----------------------------------------------------------------------- %
output = solve_mpc(initial, params, cost, constraints, functions, ...
                   0, target_ref);
% output = solve_slack_mpc(initial, params, cost, constraints, functions, ...
%                    0, target_ref);


% ----------------------------------------------------------------------- %
% Extract results
% ----------------------------------------------------------------------- %
figure(1);
clf;
figure(2);
clf;
figure(3);
clf;

figure(1);
hold on;
grid on;
title('MPC input');
stairs(output.time, output.control)

figure(2);
sgtitle('Physical state trajectory');
subplot(2,1,1);
hold on;
grid on;
plot(output.time, output.state(1:2:end), 'k-');
ylabel('x-position (m)');
subplot(2,1,2);
hold on;
grid on;
plot(output.time, output.state(2:2:end), 'k-');
ylabel('velocity (m)');
xlabel('Time (s)');

figure(3);
hold on;
title('Output');
grid on;
plot(output.time, output.Z);
ylabel('x-position (m)');
xlabel('Time (s)');