% ======================================================================= %
%                           Constrained Ball MPC                          %
% ======================================================================= %
% This script provides an example for running the MPC solver with a simple
% 1D ball. A constant reference is defined by the "constant_reference.m"
% function, and only position is tracked.

clearvars;
deg     = pi/180;
rad     = 180/pi;

% ----------------------------------------------------------------------- %
% Define MPC Parameters
% ----------------------------------------------------------------------- %
params.timestep             = 0.25;
params.horizon              = 50;
params.sim_time             = 30;

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
cost_weightings.output      = eye(1);
cost_weightings.control     = eye(1);

penalty_method              = 'linear';
penalty_weight              = 1e3;          % For slack variables
constraint_type             = 'none'; % Options are none, soft, hard, mixed

constraints.hard.rate       = [-1.5, 1.5];
constraints.hard.input      = [-1.5, 1.5];
constraints.hard.output     = [-2, 8];

constraints.soft.rate       = [-1, 1];
constraints.soft.input      = [-1.2, 1.2];
constraints.soft.output     = [0, 7];

constraints.type            = constraint_type;

% ----------------------------------------------------------------------- %
% Construct MPC Input
% ----------------------------------------------------------------------- %
mpc_input.control_model     = control_model;
mpc_input.cost              = cost_weightings;
mpc_input.constraints       = constraints;
mpc_input.penalty_method    = penalty_method;
mpc_input.penalty_weight    = penalty_weight;
mpc_input.params            = params;
mpc_input.reference         = target_ref;
mpc_input.initial           = initial;

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

