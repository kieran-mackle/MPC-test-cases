% ======================================================================= %
%                           Constrained Ball MPC                          %
% ======================================================================= %
% Run MPC on simple ball using GPOPS velocity profile output as a tracking
% reference.
% State variables: position, velocity, force rate (dynamics function shared
% with GPOPS).
% Reference function: get_reference.m

clearvars -except gpops_output
deg     = pi/180;
rad     = 180/pi;


% ----------------------------------------------------------------------- %
% Extract GPOPS output reference
% ----------------------------------------------------------------------- %
t = gpops_output.time;
x = gpops_output.state(:,1);
p = polyfit(t,x,5);
dp = polyder(p);


% ----------------------------------------------------------------------- %
% Define MPC Parameters
% ----------------------------------------------------------------------- %
params.timestep             = 0.1;
params.horizon              = 100;
params.sim_time             = 20;

% ----------------------------------------------------------------------- %
% Define function handles
% ----------------------------------------------------------------------- %
functions.dynamics          = @ball_dynamics;
functions.output            = @ball_output;

% ----------------------------------------------------------------------- %
% Define initial state
% ----------------------------------------------------------------------- %
x0                          = [0, 0, 0];
u0                          = [0.0]';

initial.state               = x0;
initial.control             = u0;

set_points                  = dp; %[5, 0]';

% Define control model
control_model.auxdata.mass  = 10;
control_model.dynamics      = @ball_dynamics;
control_model.output        = @ball_output;

% Define plant model
plant_model                 = control_model;
plant_model.auxdata.mass    = 10;

% ----------------------------------------------------------------------- %
% Define cost and constraint matrices
% ----------------------------------------------------------------------- %
cost_weightings.output      = [1, 0;
                               0, 0];
cost_weightings.control     = eye(1);

penalty_method              = 'linear';
penalty_weight              = 1e3;          % For slack variables
constraint_type             = 'none'; % Options are none, soft, hard, mixed

constraints.hard.rate       = [-1.5, 1.5];      % Force acceleration constraint
constraints.hard.input      = [-1.5, 1.5];      % Force rate constraints
constraints.hard.output     = [-2, 8;           % Position constraints
                               -5, 5];          % Velocity constraints

constraints.soft.rate       = [-1, 1];
constraints.soft.input      = [-1.2, 1.2];
constraints.soft.output     = [-1, 7;
                               -4, 4];

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
mpc_input.set_points        = set_points;
mpc_input.initial           = initial;
mpc_input.t                 = 0;

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
input.reference_function    = @get_reference;

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

% INPUTS
% --------------------------
figure(1);
set(gcf,'color','w');
sgtitle('MPC input');
subplot(2,1,1);
hold on;
grid on;
ylabel('Input force (N)');
stairs(output.time, output.state(3,:), 'k-');
stairs(gpops_output.time, gpops_output.state(:,3), 'b-');
legend('MPC', 'GPOPS');

subplot(2,1,2);
hold on;
grid on;
stairs(output.time, output.control, 'k-');
stairs(gpops_output.time, gpops_output.control, 'b-');
legend('MPC', 'GPOPS');
ylabel('Force rate (N/s)');
xlabel('Time (s)');


% STATES
% --------------------------
figure(2);
set(gcf,'color','w');
sgtitle('Physical state trajectory');
subplot(2,1,1);
hold on;
grid on;
plot(output.time, output.state(1,:), 'k-');
ylabel('x-position (m)');
subplot(2,1,2);
hold on;
grid on;
plot(output.time, output.state(2,:), 'k-');
ylabel('velocity (m)');
xlabel('Time (s)');

% OUTPUTS
% --------------------------
figure(3);
set(gcf,'color','w');
subplot(2,1,1);
hold on; grid on;
title('Position output');
plot(output.time, output.Z(1,:), 'k-');
stairs(output.time, output.ref(1,:), 'r--');
plot(t, x, 'b-');
legend('MPC response', 'Reference', 'GPOPS manouevre');
ylabel('x-position (m)');

subplot(2,1,2);
hold on;
title('Velocity output');
grid on;
plot(output.time, output.Z(2,:), 'k-');
stairs(output.time, output.ref(2,:), 'b--');
% stairs(output.time, polyval(dp, output.time));
xlabel('Time (s)');
legend('MPC response', 'Reference');


