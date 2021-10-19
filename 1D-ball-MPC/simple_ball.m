% ======================================================================= %
%                           Constrained Ball MPC                          %
% ======================================================================= %
% This script provides an example for running the MPC solver with a simple
% 1D ball. A constant reference is defined by the "constant_reference.m"
% function, and only position is tracked.

clearvars;
deg     = pi/180;
rad     = 180/pi;
tic;
% ----------------------------------------------------------------------- %
% Define MPC Parameters
% ----------------------------------------------------------------------- %
params.timestep             = 0.25;
params.horizon              = 75;
params.sim_time             = 30;
convex_solver               = 'gurobi';

% ----------------------------------------------------------------------- %
% Define function handles
% ----------------------------------------------------------------------- %
functions.dynamics          = @ball_dynamics;
functions.output            = @ball_output;

% ----------------------------------------------------------------------- %
% Define initial state
% ----------------------------------------------------------------------- %
x0                          = [0, 0];
u0                          = [0.0]';

initial.state               = x0;
initial.control             = u0;

target_ref                  = [5]';

% Define control model
control_model.auxdata.mass  = 10;
control_model.dynamics      = @ball_dynamics;
control_model.output        = @ball_output;

% Define plant model
plant_model                 = control_model;
% plant_model.auxdata.mass    = 6;

% ----------------------------------------------------------------------- %
% Define cost and constraint matrices
% ----------------------------------------------------------------------- %
Q       = 100*eye(1); % Output cost
R       = 1*eye(1); % Control cost

constraints.hard.rate       = [0, 0];
constraints.hard.input      = [0, 0];
constraints.hard.output     = [0, 0];

constraints.weights.hard_rate = [0,0];
constraints.weights.hard_input = [0,0];
constraints.weights.hard_output = [0,0];

% ----------------------------------------------------------------------- %
% Define scaling matrices
% ----------------------------------------------------------------------- %
u_nom = [5];    % Nominal control input value
p_nom = [5];    % Nominal output value

% ----------------------------------------------------------------------- %
% Scale cost matrices
% ----------------------------------------------------------------------- %
cost_weightings.control    = u_nom' * R * u_nom;
cost_weightings.output     = p_nom' * Q * p_nom;

% ----------------------------------------------------------------------- %
% Construct MPC Input
% ----------------------------------------------------------------------- %
mpc_input.control_model     = control_model;
mpc_input.cost              = cost_weightings;
mpc_input.constraints       = constraints;
% mpc_input.penalty_method    = penalty_method;
% mpc_input.penalty_weight    = penalty_weight;
mpc_input.params            = params;
mpc_input.reference         = target_ref;
mpc_input.initial           = initial;
mpc_input.solver            = convex_solver;

% ----------------------------------------------------------------------- %
% Construct Forward Simulation Input
% ----------------------------------------------------------------------- %
dyn_input.phase.state       = initial.state;
sim_input.dynamics_input    = dyn_input;
sim_input.plant_model       = plant_model;
sim_input.Ts                = params.timestep;

% ----------------------------------------------------------------------- %
% Construct Input
% ----------------------------------------------------------------------- %
input.mpc_input             = mpc_input;
input.sim_input             = sim_input;
input.reference_function    = @constant_reference;

% ----------------------------------------------------------------------- %
% Solve MPC QP Problem
% ----------------------------------------------------------------------- %
output                      = mpc_control(input);
% output                      = scaled_mpc_control(input);

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
set(gcf,'color','w');
hold on;
grid on;
title('MPC input');
ylabel('Input force (N)');
xlabel('Time (s)');
stairs(output.time, output.control, 'k-');
% patch('vertices', [0, constraints.soft.input(1); 
%                    0, constraints.hard.input(1); 
%                    params.sim_time, constraints.hard.input(1); 
%                    params.sim_time, constraints.soft.input(1)], ...
%     'faces', [1, 2, 3, 4], ...
%     'FaceColor', 'r', ...
%     'FaceAlpha', 0.3, ...
%     'LineStyle', 'none');
% plot([0,params.sim_time], [constraints.hard.input(1), ...
%      constraints.hard.input(1)], 'r-', 'LineWidth', 2);
% patch('vertices', [0, constraints.soft.input(2); 
%                    0, constraints.hard.input(2); 
%                    params.sim_time, constraints.hard.input(2); 
%                    params.sim_time, constraints.soft.input(2)], ...
%     'faces', [1, 2, 3, 4], ...
%     'FaceColor', 'r', ...
%     'FaceAlpha', 0.3, ...
%     'LineStyle', 'none');
% plot([0,params.sim_time], [constraints.hard.input(2), ...
%      constraints.hard.input(2)], 'r-', 'LineWidth', 2);

figure(2);
set(gcf,'color','w');
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
set(gcf,'color','w');
hold on;
title('Output');
grid on;
plot(output.time, output.Z, 'k-');
% plot([0,params.sim_time], [target_ref, target_ref], 'b--');
% patch('vertices', [0, constraints.soft.output(1); 
%                    0, constraints.hard.output(1); 
%                    params.sim_time, constraints.hard.output(1); 
%                    params.sim_time, constraints.soft.output(1)], ...
%     'faces', [1, 2, 3, 4], ...
%     'FaceColor', 'r', ...
%     'FaceAlpha', 0.3, ...
%     'LineStyle', 'none');
% plot([0,params.sim_time], [constraints.hard.output(1), ...
%      constraints.hard.output(1)], 'r-', 'LineWidth', 2);
% patch('vertices', [0, constraints.soft.output(2); 
%                    0, constraints.hard.output(2); 
%                    params.sim_time, constraints.hard.output(2); 
%                    params.sim_time, constraints.soft.output(2)], ...
%     'faces', [1, 2, 3, 4], ...
%     'FaceColor', 'r', ...
%     'FaceAlpha', 0.3, ...
%     'LineStyle', 'none');
% plot([0,params.sim_time], [constraints.hard.output(2), ...
%      constraints.hard.output(2)], 'r-', 'LineWidth', 2);
ylabel('x-position (m)');
xlabel('Time (s)');

toc
